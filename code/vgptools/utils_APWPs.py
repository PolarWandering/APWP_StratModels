from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from pmagpy import pmag, ipmag

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.geodesic import Geodesic
from shapely.geometry import Polygon

from vgptools.utils import spherical2cartesian, shape, eigen_decomposition, cartesian2spherical, GCD_cartesian, get_angle, PD

from joblib import Parallel, delayed


def circular_mean(longitudes: Iterable[float]) -> float:
    """Return the circular mean of longitudes in degrees on ``[0, 360)``."""

    lon_rad = np.deg2rad(np.asarray(longitudes, dtype=float))
    mean_lon = np.arctan2(np.nanmean(np.sin(lon_rad)), np.nanmean(np.cos(lon_rad)))
    return float(np.rad2deg(mean_lon) % 360.0)


def align_longitudes(longitudes: Iterable[float], center: float | None = None) -> np.ndarray:
    """Unwrap longitudes onto the 360-degree branch centered on ``center``."""

    lon = np.asarray(longitudes, dtype=float)
    if center is None:
        center = circular_mean(lon)
    return center + ((lon - center + 180.0) % 360.0 - 180.0)


def cartesian_mean_path(
    paths: pd.DataFrame,
    age_col: str = "age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
) -> pd.DataFrame:
    """Summarize an ensemble of paleopoles at each age on the unit sphere."""

    rows: list[dict[str, float]] = []
    for age, group in paths.groupby(age_col):
        xyz = np.array(
            [
                spherical2cartesian(
                    [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
                )
                for _, row in group.iterrows()
            ]
        )
        mean_xyz = xyz.mean(axis=0)
        mean_xyz /= np.linalg.norm(mean_xyz)
        latitude = np.rad2deg(np.arcsin(mean_xyz[2]))
        longitude = np.rad2deg(np.arctan2(mean_xyz[1], mean_xyz[0])) % 360.0
        aligned_lon = align_longitudes(group[longitude_col].to_numpy(), longitude)
        rows.append(
            {
                "age": float(age),
                "plon": longitude,
                "plat": latitude,
                "plon_q025": np.nanquantile(aligned_lon, 0.025),
                "plon_q50": np.nanquantile(aligned_lon, 0.5),
                "plon_q975": np.nanquantile(aligned_lon, 0.975),
                "plat_q025": np.nanquantile(group[latitude_col], 0.025),
                "plat_q50": np.nanquantile(group[latitude_col], 0.5),
                "plat_q975": np.nanquantile(group[latitude_col], 0.975),
            }
        )
    return pd.DataFrame(rows).sort_values("age").reset_index(drop=True)


def summarize_quantiles(
    values: pd.DataFrame,
    age_col: str = "age",
    value_col: str = "APW_rate",
) -> pd.DataFrame:
    """Calculate median and empirical 95 percent limits by age."""

    return (
        values.groupby(age_col)[value_col]
        .quantile([0.025, 0.5, 0.975])
        .unstack()
        .rename(columns={0.025: "q025", 0.5: "q50", 0.975: "q975"})
        .reset_index()
        .rename(columns={age_col: "age"})
        .sort_values("age")
        .reset_index(drop=True)
    )


def calculate_coherent_path_rates(
    paths: pd.DataFrame,
    path_col: str = "run",
    age_col: str = "age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
    output_col: str = "APW_rate",
) -> pd.DataFrame:
    """Calculate APW rates while preserving each sampled path trajectory.

    Adjacent poles are always taken from the same bootstrap or posterior
    realization. This is the rate estimator used for the manuscript figures.
    """

    frames: list[pd.DataFrame] = []
    for _, group in paths.sort_values([path_col, age_col]).groupby(path_col):
        group = group.copy()
        vectors = [
            spherical2cartesian(
                [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
            )
            for _, row in group.iterrows()
        ]
        distances = [np.nan]
        for index in range(1, len(vectors)):
            distances.append(np.degrees(GCD_cartesian(vectors[index - 1], vectors[index])))
        group[output_col] = np.asarray(distances) / group[age_col].diff().to_numpy()
        frames.append(group)
    return pd.concat(frames, ignore_index=True)


def gallo_effective_age_resample(
    ensemble: pd.DataFrame,
    ages: Sequence[float] | np.ndarray,
    n_paths: int,
    random_seed: int | None = None,
    effective_age_col: str = "effective_age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
) -> pd.DataFrame:
    """Reproduce the Gallo et al. (2023) pooled effective-age rate estimator.

    For every requested rounded effective age, one pole is sampled from all
    ensemble rows assigned to that age. Sampling is independent at neighboring
    ages, so synthetic paths generally combine poles from different original
    bootstrap runs.
    """

    requested_ages = np.asarray(ages, dtype=float)
    if n_paths <= 0:
        raise ValueError("n_paths must be positive")
    if requested_ages.size < 2:
        raise ValueError("At least two effective ages are required")
    if np.any(np.diff(requested_ages) <= 0):
        raise ValueError("Effective ages must be unique and strictly increasing")

    pools: dict[float, pd.DataFrame] = {}
    for age in requested_ages:
        pool = ensemble.loc[np.isclose(ensemble[effective_age_col], age)].reset_index(drop=True)
        if pool.empty:
            raise ValueError(f"No ensemble rows have {effective_age_col}={age:g}")
        pools[float(age)] = pool

    rng = np.random.default_rng(random_seed)
    frames: list[pd.DataFrame] = []
    for synthetic_run in range(n_paths):
        rows = []
        for age in requested_ages:
            pool = pools[float(age)]
            rows.append(pool.iloc[int(rng.integers(len(pool)))].copy())
        path = pd.DataFrame(rows).reset_index(drop=True)
        if "run" in path:
            path = path.rename(columns={"run": "source_run"})
        path[effective_age_col] = requested_ages
        path["run"] = synthetic_run

        vectors = [
            spherical2cartesian(
                [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
            )
            for _, row in path.iterrows()
        ]
        distances = [np.nan]
        for index in range(1, len(vectors)):
            distances.append(np.degrees(GCD_cartesian(vectors[index - 1], vectors[index])))
        path["GCD"] = distances
        path["APW_rate"] = path["GCD"] / path[effective_age_col].diff()
        frames.append(path)

    return pd.concat(frames, ignore_index=True)


def running_mean_APWP(data, plon_label = 'plon', plat_label='plat', age_label = 'age',
                                 window_length=20, time_step=1, max_age=65, min_age=0):
    """
    Returns a data frame with a running mean (Moving average) APWP..
    
    Parameters: 
    - Data Frame with the folowing columns that need to be set [['vgp_lat'],['vgp_lon'],['age']]
    - time-step of the moving average (in Ma)
    - window size of the moving average (in Ma)
    
    * it also calculates descriptive parameters for the underlying distribution of VGPs within each window, e.g. : 
    - Number of VGPs,
    - Angular dispersion
    - Concentration
    - Shape (expressed as foliation, lineation, coplanarity and collinearity)
     
    * the APWP can also be described with some parameters, e.g.:
    - Apparent Polar wander rate (in degrees per Ma) between each time interval

    """
    
    mean_pole_ages = np.arange(min_age, max_age + time_step, time_step)
    
    running_means = pd.DataFrame(columns=['age','N','n_studies','k','A95','csd','plon','plat', 'foliation','lineation','collinearity','coplanarity','elong_dir',
                                         'effective_age','effective_age_std','effective_age_median','distance2age'])
    
    for age in mean_pole_ages:
        window_min = age - (window_length / 2.)
        window_max = age + (window_length / 2.)
        poles = data.loc[(data[age_label] >= window_min) & (data[age_label] <= window_max)]
        
        if poles.empty: continue
        
        number_studies = len(poles['Study'].unique())
        mean = ipmag.fisher_mean(dec=poles[plon_label].tolist(), inc=poles[plat_label].tolist())
        
        effective_age_mean = np.round(poles[age_label].to_numpy().mean()) #
        effective_age_sd = poles[age_label].to_numpy().std()
        distance2age = np.round(np.random.normal(effective_age_mean, effective_age_sd) - age)
        
        effective_age_median = np.round(np.median(poles[age_label].to_numpy()))
        
        ArrayXYZ = np.array([spherical2cartesian([np.radians(i[plat_label]), np.radians(i[plon_label])]) for _,i in poles.iterrows()])        
        if len(ArrayXYZ) > 3:
            shapes = shape(ArrayXYZ)
            PrinComp=PD(ArrayXYZ)
            eVal, eVec = eigen_decomposition(ArrayXYZ)
            elong_dir = np.degrees(cartesian2spherical(eVec[:,1]))[1] # from T&K2004 (declination od the intermediate Evec)
            # mean['inc']=np.degrees(cartesian2spherical(PrinComp))[0]
            # mean['dec']=np.degrees(cartesian2spherical(PrinComp))[1]
        else:
            shapes = [np.nan,np.nan,np.nan,np.nan]
        
        if len(poles)>2: #ensures that dict isn't empty
            running_means.loc[age] = [age, mean['n'], number_studies, mean['k'],mean['alpha95'], mean['csd'], mean['dec'], mean['inc'], 
                                      shapes[0], shapes[1], shapes[2], shapes[3], elong_dir,
                                     effective_age_mean, effective_age_sd, effective_age_median, distance2age]
    # Set longitudes in [-180, 180]
    running_means['plon'] = running_means.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)   
    
    # The following block calculates rate of polar wander (degrees per million years) 
    running_means['PPcartesian'] = running_means.apply(lambda row: spherical2cartesian([np.radians(row['plat']),np.radians(row['plon'])]), axis = 1)
    running_means['PP_prev'] = running_means['PPcartesian'].shift(periods = 1)
    running_means['PP_next'] =  running_means['PPcartesian'].shift(periods = -1)
    running_means['GCD'] = running_means.apply(lambda row: np.degrees(GCD_cartesian(row['PP_prev'], row['PPcartesian'])), axis = 1)
    running_means['APW_rate'] = running_means['GCD']/running_means['age'].diff()

    running_means['APW_rate_eff_age'] = running_means['GCD']/running_means['effective_age'].diff().replace(0, np.nan, inplace=False)
    
    # Calculate a 'kink' angle for each position of the path
    running_means['angle'] = running_means.apply(lambda row: get_angle(row['PP_prev'], row['PPcartesian'], row['PP_next']), axis = 1)

    running_means = running_means.drop(['PPcartesian', 'PP_prev', 'PP_next'], axis=1)      
    running_means.reset_index(drop=1, inplace=True)
    
    #set the present day field for the present
    running_means['plat'] = np.where(running_means['age']==0, -90, running_means['plat'])
    running_means['plon'] = np.where(running_means['age']==0, 0, running_means['plon'])
    
    return running_means

def weighted_moving_average_APWP(data, plon_label = 'plon', plat_label='plat', age_label = 'age', study_label='study', window_length=20, time_step=1, max_age=1115, min_age=1075):
    """
    Returns a data frame with a running mean (Moving average) APWP..
    
    Parameters: 
    - Data Frame with the folowing columns that need to be set [['vgp_lat'],['vgp_lon'],['age']]
    - time-step of the moving average (in Ma)
    - window size of the moving average (in Ma)
    
    * it also calculates descriptive parameters for the underlying distribution of VGPs within each window, e.g. : 
    - Number of VGPs,
    - Angular dispersion
    - Concentration
    - Shape (expressed as foliation, lineation, coplanarity and collinearity)
     
    * the APWP can also be described with some parameters, e.g.:
    - Apparent Polar wander rate (in degrees per Ma) between each time interval

    """
    
    mean_pole_ages = np.arange(min_age, max_age + time_step, time_step)
    
    running_means = pd.DataFrame(columns=['age','N','n_studies','k','A95','csd','plon','plat', 'foliation','lineation','collinearity','coplanarity','elong_dir',
                                         'effective_age','effective_age_std','effective_age_median','distance2age'])
    
    for age in mean_pole_ages:
        window_min = age - (window_length / 2.)
        window_max = age + (window_length / 2.)
        poles_ = data.loc[(data[age_label] >= window_min) & (data[age_label] <= window_max)]
        
        weights = [(1-(np.abs( row.age - age ) / ((window_max - window_min) / 2))) for i, row in poles_.iterrows()]
        poles = poles_.sample(n = len(poles_), weights = weights, replace = True)
        
        if poles.empty: continue
        
        number_studies = len(poles[study_label].unique())
        mean = ipmag.fisher_mean(dec=poles[plon_label].tolist(), inc=poles[plat_label].tolist())
        # print(mean)
        effective_age_mean = np.round(poles[age_label].to_numpy().mean()) #
        effective_age_sd = poles[age_label].to_numpy().std()
        distance2age = np.round(np.random.normal(effective_age_mean, effective_age_sd) - age)
        
        effective_age_median = np.round(np.median(poles[age_label].to_numpy()))
        
        ArrayXYZ = np.array([spherical2cartesian([np.radians(i[plat_label]), np.radians(i[plon_label])]) for _,i in poles.iterrows()])        
        elong_dir = np.nan
        if len(ArrayXYZ) > 3:
            shapes = shape(ArrayXYZ)
            PrinComp=PD(ArrayXYZ)
            eVal, eVec = eigen_decomposition(ArrayXYZ)
            elong_dir = np.degrees(cartesian2spherical(eVec[:,1]))[1] # from T&K2004 (declination of the intermediate Evec)
            

        else:
            shapes = [np.nan,np.nan,np.nan,np.nan]
        
        if len(poles)>2: #ensures that dict isn't empty
            running_means.loc[age] = [age, mean['n'], number_studies, mean['k'],mean['alpha95'], mean['csd'], mean['dec'], mean['inc'], 
                                      shapes[0], shapes[1], shapes[2], shapes[3], elong_dir,
                                     effective_age_mean, effective_age_sd, effective_age_median, distance2age]
    # Set longitudes in [-180, 180]
    # running_means['plon'] = running_means.apply(lambda row: row.plon - 360 if row.plon > 180 else row.plon, axis =1)   
    
    # The following block calculates rate of polar wander (degrees per million years) 
    running_means['PPcartesian'] = running_means.apply(lambda row: spherical2cartesian([np.radians(row['plat']),np.radians(row['plon'])]), axis = 1)
    running_means['PP_prev'] = running_means['PPcartesian'].shift(periods = 1)
    running_means['PP_next'] =  running_means['PPcartesian'].shift(periods = -1)
    running_means['GCD'] = running_means.apply(lambda row: np.degrees(GCD_cartesian(row['PP_prev'], row['PPcartesian'])), axis = 1)
    running_means['APW_rate'] = running_means['GCD']/running_means['age'].diff()

    running_means['APW_rate_eff_age'] = running_means['GCD']/running_means['effective_age'].diff().replace(0, np.nan, inplace=False)
    
    # Calculate a 'kink' angle for each position of the path
    running_means['angle'] = running_means.apply(lambda row: get_angle(row['PP_prev'], row['PPcartesian'], row['PP_next']), axis = 1)

    running_means = running_means.drop(['PPcartesian', 'PP_prev', 'PP_next'], axis=1)      
    running_means.reset_index(drop=1, inplace=True)
    
    return running_means



def get_pseudo_vgps(df):  
    '''
    takes a DF with paleomagnetic poles and respective statistics, it draws N randomly generated VGPs
    following the pole location and kappa concentration parameter. In the present formulation we follow
    a very conservative apporach for the assignaiton of ages to each VGP, it is taken at random between
    the lower and upper bounds of the distribution of reported VGPs.
    Note: column labels are presently hard-coded into this, if relevant.
    '''
    
    data = {'Study': [], 'Plat': [], 'Plon': [], 'mean_age': []}

    for index, row in df.iterrows():
        directions_temp = ipmag.fishrot(k = row.K, n = row.N, dec = row.Plon, inc = row.Plat, di_block = False)
        
        vgp_lon_bst = directions_temp[0]
        vgp_lat_bst = directions_temp[1]
    
        if row.uncer_dist == 'uniform':
            ages = [np.random.randint(np.floor(row.min_age),np.ceil(row.max_age)) for _ in range(row.N)]    
        elif row.uncer_dist == 'normal':
            ages = [np.random.normal(row.mean_age,(row.max_age - row.min_age) / 2) for _ in range(row.N)]
            
        studies = [row.Study for _ in range(row.N)]
 
        data['Study'] += studies
        data['Plat'] += vgp_lat_bst
        data['Plon'] += vgp_lon_bst
        data['mean_age'] += ages
    
    pseudo_vgps = pd.DataFrame(data)

    return pseudo_vgps



def get_vgps_sampling_from_direction(df, study_label= 'pole_name',
                                     site_label='site',
                                     height_label='height',
                                     slat_label='lat', slon_label='lon', 
                                     n_label='dir_n_samples', 
                                     dec_label='dir_dec', inc_label='dir_inc', k_label='dir_k',
                                     f_label='dir_f', f_sigma_label='dir_f_sigma',
                                     mean_age_lab='age', min_age_lab='age_low', max_age_lab='age_high', 
                                     age_uncertainty_label='age_uncertainty', polarity_label='dir_polarity'):
    
    '''
    Input:
    
    DF with the following site-level information: 
    - study, site coordinates, mean direction, concentration parameter, mean age and error distribution. 
    
    Steps:
    1. generate a pseudo-sample from the original dataset using nonparametric random sampling 
    with replacement (bootstrap sample).
    2. For each row (site-level entry) in the bootstrap sample draws a random direction following 
    the kappa concentration parameter and mean direction. Note here our assumption is that site-level uncertainties 
    are associated with random errors generated during field sample collection and lab sample handling and during 
    experiements such that they can be treated as Fisher distributed around the mean directions. Note that this 
    assumed error circle in directional space would not map to a circularly symmetric error ellipse in pole space.
    3. Assign an age to the resampled directional data from the corresponding age assignments 
    and associated type of uncertainty (uniform or gaussian).
    3. For each row, calculates the corresponding VGP. 
    
    Note: input directions must be all in the same mode (sensu PmagPy) so that we get coherent vgps. 
    Given that we are used to looking at the Keweenawan Track in the northern hemisphere, we will work with normal inclinations   
    
    Output:
    
    - A DataFrame with the same size than the original dataset, with randomized parameters.  
    '''    

    study, site, height, age_bst, decs, incs, slat, slon, spolarity, indexes = [], [], [], [], [], [], [], [], [], [] # parameters of the pseudo-sample to be filled

    for index, row in df.iterrows():
        # in case the record is sedimentary (Freda or Nonesuch) each sample is a site
        if row['geologic_classes'] == 'Sedimentary':
            # instead of resampling the site direction with a madeup kappa, we consider the uncertainty in f factor to resample the inclination-corrected site directions
            f_nominal = row[f_label]
            f_sigma = row[f_sigma_label]
            resample_f = np.random.normal(f_nominal, f_sigma)
            decs.append(row[dec_label])
            incs.append(pmag.unsquish(row[inc_label], resample_f))
        else:
            # else we are dealing with igneous site level records
            # we first generate one random direction from the original entry.
            kappa = row[k_label] 

            # make sure to generate at least three sample directions to take a site mean    
            n = 3 if (np.isnan(row[n_label])) or (row[n_label]<3) else int(row[n_label])

            # draw random directions from the Fisher distribution
            directions_temp = ipmag.fishrot(k = kappa, n = n, dec = row[dec_label], inc = row[inc_label], di_block = False)    
            # take a site mean    
            site_mean = ipmag.fisher_mean(dec = directions_temp[0], inc = directions_temp[1])        
            
            decs.append(site_mean['dec'])
            incs.append(site_mean['inc'])

        site.append(row[site_label])
        height.append(row[height_label])
        slat.append(row[slat_label])
        slon.append(row[slon_label])
        indexes.append(index)
        study.append(row[study_label])
        spolarity.append(row[polarity_label])

        # Assessing the age uncertianty distribution (uniform or normal)
        if row[age_uncertainty_label] == 'uniform':
            age_bst.append(np.random.uniform(np.floor(row[min_age_lab]),np.ceil(row[max_age_lab])))
        else:            
            age_bst.append(np.random.normal(row[mean_age_lab], (row[max_age_lab] - row[mean_age_lab]) / 2)) 
    
    dictionary = {
                  'study': study,
                  'site': site,
                  'height': height,
                  'age': age_bst,
                  'dec': decs,    
                  'inc': incs,
                  'slat': slat,
                  'slon': slon,
                  'spolarity': spolarity
                  }    
    new_df = pd.DataFrame(dictionary)        
    
    # calculate the corresponding VGPs
    new_VGPs = pmag.dia_vgp(new_df['dec'], new_df['inc'], 1, new_df['slat'], new_df['slon'])
    new_df['plon'] = new_VGPs[0]
    
    new_df['plat'] = new_VGPs[1]
    
    new_df['plon'] = np.where(new_df['spolarity']== 'n', new_df['plon'], (new_df['plon']+180) %360)
    new_df['plat'] = np.where(new_df['spolarity']== 'n', new_df['plat'], -new_df['plat'])

    new_df.index = indexes

    return new_df

def ensemble_RMs(df_vgps_original, n_sims = 100, 
                 study_label= 'pole_name', slat_label='lat', slon_label='lon', 
                 dec_label='dir_dec', inc_label='dir_inc', k_label='dir_k',
                 mean_age_lab='age', min_age_lab='age_low', max_age_lab='age_high',
                 plon_label = 'plon', plat_label='plat', age_label = 'age',
                 window_length=20, time_step=1, max_age=1115, min_age=1075):
    
    '''
    This function gets the pooled directions and return the ensemble Moving Averages following the workflow of Gallo et al. (2023) 
    '''
    
    running_means_global = pd.DataFrame(columns=['run','N','k','A95','csd','foliation','lineation','collinearity','coplanarity'])
    
    # make an ensemble of pseudo-dataframes for post-processing age distribution
    ensemble_pseudo_df = []
    
    for i in range(n_sims):
    
        # Generate a pseudo-sample of the original dataset in which every entry is a pseudo-sample taken for the error PDF
        pseudo_df = get_vgps_sampling_from_direction(df_vgps_original, study_label= study_label,
                                             slat_label=slat_label, slon_label=slon_label, 
                                             dec_label=dec_label, inc_label=inc_label, k_label=k_label,
                                             mean_age_lab=mean_age_lab, min_age_lab=min_age_lab, max_age_lab=max_age_lab)
        pseudo_df['run'] = float(i)
        ensemble_pseudo_df.append(pseudo_df)
        # Construct a Moving Average on the former data-set
        RM = weighted_moving_average_APWP(pseudo_df, plon_label = plon_label, plat_label=plat_label, age_label = age_label, 
                            window_length=window_length, time_step=time_step, max_age=max_age, min_age=min_age)
 
        RM['run'] = float(i)
        running_means_global = pd.concat([running_means_global, RM], ignore_index=True)

    # running_means_global['plon'] = running_means_global['plon'].where(running_means_global['plon'] <= 180, running_means_global['plon'] - 360)
    ensemble_pseudo_df = pd.concat(ensemble_pseudo_df).reset_index(drop=True)
    return running_means_global, ensemble_pseudo_df



def _one_ensemble_run(i, df_vgps_original,
                      study_label, slat_label, slon_label,
                      dec_label, inc_label, k_label,
                      mean_age_lab, min_age_lab, max_age_lab,
                      plon_label, plat_label, age_label,
                      window_length, time_step, max_age, min_age):

    pseudo_df = get_vgps_sampling_from_direction(
        df_vgps_original,
        study_label=study_label,
        slat_label=slat_label,
        slon_label=slon_label,
        dec_label=dec_label,
        inc_label=inc_label,
        k_label=k_label,
        mean_age_lab=mean_age_lab,
        min_age_lab=min_age_lab,
        max_age_lab=max_age_lab
    )

    pseudo_df['run'] = float(i)

    RM = weighted_moving_average_APWP(
        pseudo_df,
        plon_label=plon_label,
        plat_label=plat_label,
        age_label=age_label,
        window_length=window_length,
        time_step=time_step,
        max_age=max_age,
        min_age=min_age
    )

    RM['run'] = float(i)

    return RM, pseudo_df

def ensemble_RMs_parallel(df_vgps_original, n_sims=100,
                 study_label='pole_name', slat_label='lat', slon_label='lon',
                 dec_label='dir_dec', inc_label='dir_inc', k_label='dir_k',
                 mean_age_lab='age', min_age_lab='age_low', max_age_lab='age_high',
                 plon_label='plon', plat_label='plat', age_label='age',
                 window_length=20, time_step=1, max_age=1115, min_age=1075,
                 n_jobs=-1, verbose=5):

    results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_one_ensemble_run)(
            i,
            df_vgps_original,
            study_label, slat_label, slon_label,
            dec_label, inc_label, k_label,
            mean_age_lab, min_age_lab, max_age_lab,
            plon_label, plat_label, age_label,
            window_length, time_step, max_age, min_age
        )
        for i in range(n_sims)
    )

    RMs, pseudo_dfs = zip(*results)

    running_means_global = pd.concat(RMs, ignore_index=True)
    ensemble_pseudo_df = pd.concat(pseudo_dfs, ignore_index=True)

    return running_means_global, ensemble_pseudo_df


def ultimate_VGP_resample(df, age_model_posterior_df, seds_f_df, 
                        study_label= 'pole_name',
                        site_label='site',
                        height_label='height',
                        height_upper_label='height_upper',
                        height_lower_label='height_lower',
                        slat_label='lat', slon_label='lon', 
                        n_label='dir_n_samples', 
                        dec_label='dir_dec', inc_label='dir_inc', k_label='dir_k', f_label = 'dir_f',
                        mean_age_lab='age', min_age_lab='age_low', max_age_lab='age_high', 
                        age_uncertainty_label='age_uncertainty', polarity_label='dir_polarity'):
    
    study, site, height, age_bst, decs, incs, slat, slon, spolarity, indexes = [], [], [], [], [], [], [], [], [], [] # parameters of the pseudo-sample to be filled

    # separate records with and without age model
    df_no_age_model = df[df['age_model'].isna()]
    df_with_age_model = df[~df['age_model'].isna()]
    
    # for those with age models, we need to load the Julia age model csv file and resample the ages based on heights
    for section in df_with_age_model['age_model'].unique():
        this_section = df_with_age_model[df_with_age_model['age_model'] == section]
        if 'uniform' in this_section['height_type'].unique():
            this_section_uniform = this_section[this_section['height_type'] == 'uniform']
            this_section_exact = this_section[this_section['height_type'] == 'exact']

            this_section_uniform['height'] = np.random.uniform(this_section_uniform[height_lower_label].min(), this_section_uniform[height_upper_label].max(), size=this_section_uniform.shape[0])
            this_section = pd.concat([this_section_uniform, this_section_exact]).reset_index(drop=True)

        this_age_model = age_model_posterior_df.loc[section, 'posterior_age_model']
        this_heights = this_age_model['height']
        this_age_dist = this_age_model.iloc[:, 1:].sample(axis=1, n=1).squeeze()
        this_age_resampled = np.interp(this_section['height'], this_heights, this_age_dist)
        df_with_age_model.loc[this_section.index, 'pseudo_age'] = this_age_resampled

    # for those without age models, we just need to resample the ages based on the age uncertainty
    for section in df_no_age_model['pole_name'].unique():
        this_section = df_no_age_model[df_no_age_model['pole_name'] == section]
        if this_section[age_uncertainty_label].unique()[0] == 'uniform':
            df_no_age_model.loc[this_section.index, 'pseudo_age'] = np.random.uniform(this_section[min_age_lab].min(), this_section[max_age_lab].max(), size=this_section.shape[0])
        else:
            df_no_age_model.loc[this_section.index, 'pseudo_age'] = np.random.normal(this_section[mean_age_lab], (this_section[max_age_lab] - this_section[mean_age_lab]) / 2)
    
    df_with_pseudo_age = pd.concat([df_no_age_model, df_with_age_model]).reset_index(drop=True)

    # now we resample f factors for sedimentary sites
    sed_sites = df_with_pseudo_age[df_with_pseudo_age['geologic_classes'] == 'Sedimentary'].copy()
    igneous_sites = df_with_pseudo_age[~(df_with_pseudo_age['geologic_classes'] == 'Sedimentary')].copy()
    sed_sites[f_label] = pd.to_numeric(sed_sites.get(f_label, np.nan), errors='coerce')
    igneous_sites[f_label] = 1.0

    # redraw f factors based on section name
    for section in sed_sites['pole_name'].unique():
        this_section = sed_sites[sed_sites['pole_name'] == section]
        f_pool = pd.to_numeric(seds_f_df.at[section, 'f_factors'].squeeze(), errors='coerce').dropna().to_numpy(dtype=float)
        sed_sites.loc[this_section.index, f_label] = np.random.choice(f_pool, size=this_section.shape[0]).astype(float)

    df_with_pseudo_age_f = pd.concat([igneous_sites, sed_sites]).reset_index(drop=True)
    
    # now we can resample the VGPs
    for index, row in df_with_pseudo_age_f.iterrows():
        # make sure to generate at least three sample directions to take a site mean    
        n = 3 if (np.isnan(row[n_label])) or (row[n_label]<3) else int(row[n_label])

        # draw random directions from the Fisher distribution
        directions_temp = ipmag.fishrot(k = row[k_label] , n = n, dec = row[dec_label], inc = row[inc_label], di_block = False)    
        # take a site mean    
        site_mean = ipmag.fisher_mean(dec = directions_temp[0], inc = directions_temp[1])        
        
        decs.append(site_mean['dec'])

        # apply f factor unsquishing for sedimentary sites
        incs.append(pmag.unsquish(site_mean['inc'], row[f_label]))

        site.append(row[site_label])
        height.append(row[height_label])
        slat.append(row[slat_label])
        slon.append(row[slon_label])
        indexes.append(index)
        study.append(row[study_label])
        spolarity.append(row[polarity_label])

        age_bst.append(row['pseudo_age'])
    
    dictionary = {
                  'study': study,
                  'site': site,
                  'height': height,
                  'age': age_bst,
                  'dec': decs,    
                  'inc': incs,
                  'slat': slat,
                  'slon': slon,
                  'spolarity': spolarity
                  }    
    new_df = pd.DataFrame(dictionary)        
    
    # calculate the corresponding VGPs
    new_VGPs = pmag.dia_vgp(new_df['dec'], new_df['inc'], 1, new_df['slat'], new_df['slon'])
    new_df['plon'] = new_VGPs[0]
    
    new_df['plat'] = new_VGPs[1]
    
    new_df['plon'] = np.where(new_df['spolarity']== 'n', new_df['plon'], (new_df['plon']+180) %360)
    new_df['plat'] = np.where(new_df['spolarity']== 'n', new_df['plat'], -new_df['plat'])

    new_df.index = indexes

    return new_df


def ultimate_ensemble_RMs(df_vgps_original, age_model_posterior_df, seds_f_df, n_sims = 100, 
                        study_label= 'pole_name',
                        site_label='site',
                        height_label='height',
                        height_upper_label='height_upper',
                        height_lower_label='height_lower',
                        slat_label='lat', slon_label='lon', 
                        n_label='dir_n_samples', 
                        dec_label='dir_dec', inc_label='dir_inc', k_label='dir_k', f_label = 'dir_f',
                        mean_age_lab='age', min_age_lab='age_low', max_age_lab='age_high', 
                        age_uncertainty_label='age_uncertainty', polarity_label='dir_polarity',
                        plon_label = 'plon', plat_label='plat', age_label = 'age',
                        window_length=20, time_step=1, max_age=1115, min_age=1075):
    
    running_means_global = pd.DataFrame(columns=['run','N','k','A95','csd','foliation','lineation','collinearity','coplanarity'])
    
    # make an ensemble of pseudo-dataframes for post-processing age distribution
    ensemble_pseudo_df = []
    
    for i in range(n_sims):
    
        # Generate a pseudo-sample of the original dataset in which every entry is a pseudo-sample taken for the error PDF
        pseudo_df = ultimate_VGP_resample(df_vgps_original, age_model_posterior_df, seds_f_df, 
                                          study_label= study_label,
                                          site_label=site_label, height_label=height_label, height_upper_label=height_upper_label, height_lower_label=height_lower_label,
                                          slat_label=slat_label, slon_label=slon_label, n_label=n_label,
                                          dec_label=dec_label, inc_label=inc_label, k_label=k_label,
                                          mean_age_lab=mean_age_lab, min_age_lab=min_age_lab, max_age_lab=max_age_lab,
                                          age_uncertainty_label=age_uncertainty_label, polarity_label=polarity_label, f_label=f_label)
        pseudo_df['run'] = float(i)
        ensemble_pseudo_df.append(pseudo_df)
        # Construct a Moving Average on the former data-set
        RM = weighted_moving_average_APWP(pseudo_df, plon_label = plon_label, plat_label=plat_label, age_label = age_label, 
                            window_length=window_length, time_step=time_step, max_age=max_age, min_age=min_age)
 
        RM['run'] = float(i)
        running_means_global = pd.concat([running_means_global, RM], ignore_index=True)

    # running_means_global['plon'] = running_means_global['plon'].where(running_means_global['plon'] <= 180, running_means_global['plon'] - 360)
    ensemble_pseudo_df = pd.concat(ensemble_pseudo_df).reset_index(drop=True)
    return running_means_global, ensemble_pseudo_df

def _one_ultimate_ensemble_run(
    i,
    df_vgps_original,
    age_model_posterior_df,
    seds_f_df,
    study_label,
    site_label,
    height_label,
    height_upper_label,
    height_lower_label,
    slat_label,
    slon_label,
    n_label,
    dec_label,
    inc_label,
    k_label,
    f_label,
    mean_age_lab,
    min_age_lab,
    max_age_lab,
    age_uncertainty_label,
    polarity_label,
    plon_label,
    plat_label,
    age_label,
    window_length,
    time_step,
    max_age,
    min_age,
):

    # 1. Resample VGPs
    pseudo_df = ultimate_VGP_resample(
        df_vgps_original,
        age_model_posterior_df,
        seds_f_df,
        study_label=study_label,
        site_label=site_label,
        height_label=height_label,
        height_upper_label=height_upper_label,
        height_lower_label=height_lower_label,
        slat_label=slat_label,
        slon_label=slon_label,
        n_label=n_label,
        dec_label=dec_label,
        inc_label=inc_label,
        k_label=k_label,
        mean_age_lab=mean_age_lab,
        min_age_lab=min_age_lab,
        max_age_lab=max_age_lab,
        age_uncertainty_label=age_uncertainty_label,
        polarity_label=polarity_label,
        f_label=f_label,
    )

    pseudo_df = pseudo_df.copy()
    pseudo_df["run"] = float(i)

    # 2. Running mean
    RM = weighted_moving_average_APWP(
        pseudo_df,
        plon_label=plon_label,
        plat_label=plat_label,
        age_label=age_label,
        window_length=window_length,
        time_step=time_step,
        max_age=max_age,
        min_age=min_age,
    )

    RM = RM.copy()
    RM["run"] = float(i)

    return RM, pseudo_df

def ultimate_ensemble_RMs_parallel(
    df_vgps_original,
    age_model_posterior_df,
    seds_f_df,
    n_sims=100,
    study_label="pole_name",
    site_label="site",
    height_label="height",
    height_upper_label="height_upper",
    height_lower_label="height_lower",
    slat_label="lat",
    slon_label="lon",
    n_label="dir_n_samples",
    dec_label="dir_dec",
    inc_label="dir_inc",
    k_label="dir_k",
    f_label="dir_f",
    mean_age_lab="age",
    min_age_lab="age_low",
    max_age_lab="age_high",
    age_uncertainty_label="age_uncertainty",
    polarity_label="dir_polarity",
    plon_label="plon",
    plat_label="plat",
    age_label="age",
    window_length=20,
    time_step=1,
    max_age=1115,
    min_age=1075,
    n_jobs=-1,
    backend="loky",
):
    results = Parallel(
        n_jobs=n_jobs,
        backend=backend,
        verbose=5,
    )(
        delayed(_one_ultimate_ensemble_run)(
            i,
            df_vgps_original,
            age_model_posterior_df,
            seds_f_df,
            study_label,
            site_label,
            height_label,
            height_upper_label,
            height_lower_label,
            slat_label,
            slon_label,
            n_label,
            dec_label,
            inc_label,
            k_label,
            f_label,
            mean_age_lab,
            min_age_lab,
            max_age_lab,
            age_uncertainty_label,
            polarity_label,
            plon_label,
            plat_label,
            age_label,
            window_length,
            time_step,
            max_age,
            min_age,
        )
        for i in range(n_sims)
    )

    # Unpack
    RM_list, pseudo_list = zip(*results)

    running_means_global = pd.concat(RM_list, ignore_index=True)
    ensemble_pseudo_df = pd.concat(pseudo_list, ignore_index=True)

    return running_means_global, ensemble_pseudo_df
