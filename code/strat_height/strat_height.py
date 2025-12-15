import folium
from shapely.geometry import Point, LineString
from geopy.distance import geodesic
from pyhigh import get_elevation
import numpy as np
import pmagpy.ipmag as ipmag


def assign_lava_bedding_to_sites(pmag_site_data, bedding_data, dip_dir_col='DIPD_TREND', dip_col='DIP_PLUNGE'):
    '''
    function for assigning the bedding data to the sites based on the mean dip direction and dip

    Parameters
    ----------
    pmag_site_data : pandas DataFrame
        site data
    bedding_data : pandas DataFrame
        bedding data
    dip_dir_col : str
        column name for dip direction in bedding_data
    dip_col : str
        column name for dip in bedding_data

    Returns
    -------
    pmag_site_data : pandas DataFrame
        site data with the mean bedding dip direction and dip assigned
    '''


    mean_dip_orientation = ipmag.fisher_mean(bedding_data[dip_dir_col].tolist(), bedding_data[dip_col].tolist())
    mean_dip_direction = mean_dip_orientation['dec']
    mean_dip = mean_dip_orientation['inc']
    print(mean_dip_orientation)
    pmag_site_data['bedding_dip_dir'] = mean_dip_direction
    pmag_site_data['bedding_dip'] = mean_dip
    return pmag_site_data

def plot_lava(folium_map, polygon_file, color):
    '''
    function for plotting the lava flow polygon on the folium map

    Parameters
    ----------
    folium_map : folium map object
        folium map object
    polygon_file : shapely Polygon object
        lava flow polygon
    color : str
        color of the polygon
    '''
    if polygon_file.geom_type == 'Polygon':
        locations = [(lat, lon%360) for lon, lat in polygon_file.exterior.coords]
        folium.Polygon(locations=locations, color="", fill=True, fill_color=color, fill_opacity=0.6).add_to(folium_map)
    else:
        for poly in polygon_file.geoms:
            locations = [(lat, lon%360) for lon, lat in poly.exterior.coords]
            folium.Polygon(locations=locations, color="", fill=True, fill_color=color, fill_opacity=0.6).add_to(folium_map)

def plot_site(folium_map, site_data, color, lat_col='lat', lon_col='lon', popupname='site'):
    '''
    function for plotting the sites on the folium map

    Parameters
    ----------
    folium_map : folium map object
        folium map object
    site_data : pandas DataFrame
        site data
    color : str
        color of the circle marker
    lat_col : str
        column name for latitude in site_data
    lon_col : str
        column name for longitude in site_data
    popupname : str
        column name for the popup text
    '''
    for idx, row in site_data.iterrows():
        folium.CircleMarker(
            location=[row[lat_col], row[lon_col]%360],
            popup=row[popupname], 
            radius=5,  # Size of the circle marker
            color="",  # Border color of the circle
            fill=True,
            fill_color=color,  # Fill color of the circle
            fill_opacity=0.7,
        ).add_to(folium_map)

def create_line_from_point(x, y, azimuth, distance=5000):
    """
    Create a line starting from the given point with a specified azimuth and distance.

    Parameters
    ----------
    x : float
        x coordinate of the starting point
    y : float
        y coordinate of the starting point
    azimuth : float
        azimuth in degrees
    distance : float
        distance in meters

    Returns
    -------
    line : shapely LineString object
        line from the starting point with the specified azimuth and distance
    """
    meter2degree = 1/111000
    distance=distance*meter2degree
    azimuth_rad = np.radians(azimuth)
    x_end = x + distance * np.sin(azimuth_rad)
    y_end = y + distance * np.cos(azimuth_rad)
    return LineString([Point(x, y), Point(x_end, y_end)])

def estimate_height_from_strat_top(folium_map, site_data, lava_polygon, datum,
                                   distance=5000, lat_col='lat', lon_col='lon', 
                                   dip_col='bedding_dip', dip_dir_col='bedding_dip_dir', 
                                   color='purple', multipoint_selection='min'):
    '''
    function for estimating the strat height of the site from the top of a unit

    Parameters
    ----------
    folium_map : folium map object
        folium map object
    site_data : pandas DataFrame
        site data
    lava_polygon : shapely Polygon object
        lava polygon
    datum : float
        datum for the calculated strat height
        need to assign to be the basis of a lava flow as provided by Green, 2011
    distance : float
        distance in meters to draw the line from the site toward the lava base
        for calculating the intersection point
    lat_col : str
        column name for latitude in site_data
    lon_col : str
        column name for longitude in site_data
    dip_col : str
        column name for dip in site_data
    dip_dir_col : str
        column name for dip direction in site_data
    color : str
        color of the circle marker
        
    '''
    for i, r in site_data.iterrows():
        r[lon_col] = r[lon_col] % 360
        r[dip_dir_col] = r[dip_dir_col] % 360
        site_elevation = get_elevation(r[lat_col], (r[lon_col] + 180) % 360 - 180)
        line = create_line_from_point(r[lon_col], r[lat_col], r[dip_dir_col]+180, distance=distance)
        # folium.PolyLine(locations=[(r[lat_col], r[lon_col]), (line.coords[-1][1], line.coords[-1][0])], color=color).add_to(folium_map)
        intersection = line.intersection(lava_polygon.exterior)
        if intersection.type == 'Point':
            folium.CircleMarker([intersection.y, intersection.x], 
                                radius=5,  # Size of the circle marker
                                color="",  # Border color of the circle
                                fill=True,
                                fill_color=color,  # Fill color of the circle
                                fill_opacity=0.7).add_to(folium_map)
        elif intersection.type == 'MultiPoint':
            # there are multiple points, we want the one with smallest longitude in the MultiPoint object
            if multipoint_selection == 'min':
                intersection = min(intersection.geoms, key=lambda point: point.x)
            elif multipoint_selection == 'max':
                intersection = max(intersection.geoms, key=lambda point: point.x)
            folium.CircleMarker([intersection.y, intersection.x], 
                                radius=5,  # Size of the circle marker
                                color="",  # Border color of the circle
                                fill=True,
                                fill_color=color,  # Fill color of the circle
                                fill_opacity=0.7).add_to(folium_map)
        print('site lat', r[lat_col], 'site lon', r[lon_col], 'site elevation', site_elevation, 'intersection lat', intersection.y, 'intersection lon', intersection.x)
        
        intersection_elevation = get_elevation(intersection.y, (intersection.x + 180) % 360 - 180)
        folium.PolyLine(locations=[(r[lat_col], r[lon_col]), (intersection.y, intersection.x)], color=color).add_to(folium_map)
        
        # now calculate the distance in meters of the line from the intersection point to the site
        distance_from_flowbase = geodesic((r[lat_col], r[lon_col]), (intersection.y, intersection.x)).meters

        # site_data.loc[i, 'distance_from_flowbase'] = distance_from_flowbase
        # write the distance from the flow base to the site in the original site data table
        if intersection_elevation < site_elevation:
            site_data.loc[i, 'relative height'] = (distance_from_flowbase + np.abs(site_elevation - intersection_elevation) / np.tan(np.radians(r[dip_col]))) * np.sin(np.radians(r[dip_col]))
            site_data.loc[i, 'height'] = datum + (distance_from_flowbase + np.abs(site_elevation - intersection_elevation) / np.tan(np.radians(r[dip_col]))) * np.sin(np.radians(r[dip_col]))
        else:
            site_data.loc[i, 'relative height'] = (distance_from_flowbase - np.abs(intersection_elevation - site_elevation) / np.tan(np.radians(r[dip_col]))) * np.sin(np.radians(r[dip_col]))
            site_data.loc[i, 'height'] = datum + (distance_from_flowbase - np.abs(intersection_elevation - site_elevation) / np.tan(np.radians(r[dip_col]))) * np.sin(np.radians(r[dip_col]))


# need to add optional arguments for either adding or subtracting from a given datum
def calc_strat(folium_map, lat, lon, dip_dir, dip, distance, polygon, color, multipoint_selection='min', datum=0):

    lon = lon % 360
    dip_dir = dip_dir % 360
    site_elevation = get_elevation(lat, (lon + 180) % 360 - 180)
    line = create_line_from_point(lon, lat, dip_dir, distance=distance)
    intersection = line.intersection(polygon.exterior)
    if intersection.type == 'Point':
        folium.CircleMarker([intersection.y, intersection.x], 
                            radius=5,  # Size of the circle marker
                            color="",  # Border color of the circle
                            fill=True,
                            fill_color=color,  # Fill color of the circle
                            fill_opacity=0.7).add_to(folium_map)
    elif intersection.type == 'MultiPoint':
        # there are multiple points, we want the one with smallest longitude in the MultiPoint object
        if multipoint_selection == 'min':
            intersection = min(intersection.geoms, key=lambda point: point.x)
        elif multipoint_selection == 'max':
            intersection = max(intersection.geoms, key=lambda point: point.x)
        folium.CircleMarker([intersection.y, intersection.x], 
                            radius=5,  # Size of the circle marker
                            color="",  # Border color of the circle
                            fill=True,
                            fill_color=color,  # Fill color of the circle
                            fill_opacity=0.7).add_to(folium_map)
    print('site lat', lat, 'site lon', lon, 'site elevation', site_elevation, 'intersection lat', intersection.y, 'intersection lon', intersection.x)
    
    intersection_elevation = get_elevation(intersection.y, (intersection.x + 180) % 360 - 180)
    folium.PolyLine(locations=[(lat, lon), (intersection.y, intersection.x)], color=color).add_to(folium_map)
    
    # now calculate the distance in meters of the line from the intersection point to the site
    distance_from_flowbase = geodesic((lat, lon), (intersection.y, intersection.x)).meters

    if intersection_elevation < site_elevation:
        relative_height = (distance_from_flowbase + np.abs(site_elevation - intersection_elevation) / np.tan(np.radians(dip))) * np.sin(np.radians(dip))
        height = datum - (distance_from_flowbase + np.abs(site_elevation - intersection_elevation) / np.tan(np.radians(dip))) * np.sin(np.radians(dip))
    else:
        relative_height = (distance_from_flowbase - np.abs(intersection_elevation - site_elevation) / np.tan(np.radians(dip))) * np.sin(np.radians(dip))
        height = datum - (distance_from_flowbase - np.abs(intersection_elevation - site_elevation) / np.tan(np.radians(dip))) * np.sin(np.radians(dip))

    return folium_map, relative_height, height, distance_from_flowbase