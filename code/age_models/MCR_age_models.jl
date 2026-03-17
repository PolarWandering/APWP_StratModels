function make_ChronAgeData(Name, Age, Age_sigma, Height, Age_Sidedness, Path; Age_Unit, Height_Unit)

    @assert length(Age) == length(Age_sigma) == length(Height) == length(Age_Sidedness)

    nSamples = length(Age)

    section = ChronAgeData(nSamples)
    section.Name          = Name
    section.Age           = Age # Measured ages
    section.Age_sigma     = Age_sigma # Measured 1-σ uncertainties
    section.Height        = Height # Depths below surface should be negative
    section.Height_sigma  = fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
    section.Age_Sidedness = Age_Sidedness # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    
    section.Path = Path # Where you want output files to be stored
    section.Age_Unit = "Ma"; # Unit of measurement for ages
    section.Height_Unit = "m"; # Unit of measurement for Height and Height_sigma

    return section
end

function run_age_model(section, resolution, bounding, nsteps; save_results=false, plot_results=true, save_plot=false)

    # Configure the stratigraphic Monte Carlo model
    config = StratAgeModelConfiguration()
    # If in doubt, you can probably leave these parameters as-is
    config.resolution = resolution # Same units as sample height. Smaller is slower!
    config.bounding = bounding # how far away do we place runaway bounds, as a fraction of total section height
    (bottom, top) = extrema(section.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = nsteps # Number of steps to run in distribution MCMC
    config.burnin = nsteps*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

    # Run the stratigraphic MCMC model
    (mdl, agedist, lldist) = StratMetropolis(section, config)

    if save_results
        # Write the results to file
        run(`mkdir -p $(section.Path)`) # Make sure that the path exists
        writedlm(section.Path*"_agedist.csv", agedist, ',') # Stationary distribution of the age-depth model
        writedlm(section.Path*"_height.csv", mdl.Height, ',') # Stratigraphic heights corresponding to each row of agedist
        writedlm(section.Path*"_age.csv", mdl.Age, ',') # Mean age of resulting model
        writedlm(section.Path*"_age_025CI.csv", mdl.Age_025CI, ',') # 2.5% confidence interval of resulting model
        writedlm(section.Path*"_age_975CI.csv", mdl.Age_975CI, ',') # 97.5% confidence interval of resulting model
        writedlm(section.Path*"_lldist.csv", lldist, ',') # Log likelihood distribution (to check for stationarity)
    end

    if plot_results
        # Plot the log likelihood to make sure we're converged (n.b burnin isn't recorded)
        hdl = plot(lldist,xlabel="Step number",ylabel="Log likelihood",label="",line=(0.85,:darkblue),framestyle=:box)
        if save_plot
            savefig(hdl,section.Path*"lldist.pdf")
        end
        display(hdl)
    end
    return config, mdl, agedist, lldist
end

# function plot_age_model(section, mdl; save_figure=false)
#     # Plot results (mean and 95% confidence interval for both model and data)
#     hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:grey), label="model")
#     plot!(hdl, mdl.Age, mdl.Height, linecolor=:grey, label="", fg_color_legend=:white) # Center line
#     t = section.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
#     any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=2*section.Age_sigma[t],label="data",seriestype=:scatter,color=:black)
#     t = section.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
#     any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=(2*section.Age_sigma[t],zeros(count(t))),label="",seriestype=:scatter,color=:grey,msc=:grey)
#     # any(t) && zip(section.Age[t], section.Age[t].+nanmean(section.Age_sigma[t])*4, section.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
#     t = section.Age_Sidedness .== -1 # Maximum ages (plot in orange)
#     any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=(zeros(count(t)),2*section.Age_sigma[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
#     # any(t) && zip(section.Age[t], section.Age[t].-nanmean(section.Age_sigma[t])*4, section.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
#     plot!(hdl, xlabel="Age ($(section.Age_Unit))", ylabel="Height ($(section.Height_Unit))", framestyle=:box)
#     # reverse the x axis 
#     plot!(hdl, xflip=true)
#     # put the legend at the bottom right
#     plot!(size=(200, 1000), legend=:bottomright, legendfontsize=9, xrotation=45)
#     if save_figure
#         savefig(hdl,section.Path*"AgeDepthModel.pdf")
#     end
#     return hdl
# end

function plot_age_model(section, mdl; save_figure=false)
    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.6,:grey), linecolor=:grey, label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:grey, label="", fg_color_legend=:white) # Center line
    t = section.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=2*section.Age_sigma[t], label="two-sided age", seriestype=:scatter, color=:black)
    t = section.Age_Sidedness .== 1 # Minimum ages (plot in blue with left-pointing marker)
    any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=(2*section.Age_sigma[t], zeros(count(t))), label="min age", seriestype=:scatter, color=:dodgerblue, msc=:dodgerblue, markershape=:ltriangle, markersize=6)
    # # any(t) && zip(section.Age[t], section.Age[t].+nanmean(section.Age_sigma[t])*4, section.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = section.Age_Sidedness .== -1 # Maximum ages (plot in orange with right-pointing marker)
    any(t) && plot!(hdl, section.Age[t], section.Height[t], xerror=(zeros(count(t)), 2*section.Age_sigma[t]), label="max age", seriestype=:scatter, color=:orange, msc=:orange, markershape=:rtriangle, markersize=6)
    # # any(t) && zip(section.Age[t], section.Age[t].-nanmean(section.Age_sigma[t])*4, section.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(hdl, xlabel="Age ($(section.Age_Unit))", ylabel="Stratigraphic height ($(section.Height_Unit))", framestyle=:box)
    # reverse the x axis 
    plot!(hdl, xflip=true)
    # set y axis tick font size 
    plot!(hdl, ytickfont=10, xtickfont=10)
    # put the legend at the bottom right
    plot!(size=(320, 1000), legend=:bottomright, legendfontsize=11, xrotation=45)
    if save_figure
        savefig(hdl,section.Path*"AgeDepthModel.pdf")
    end
    return hdl
end

# function plot_accumulation_model(section, config, mdl, agedist; binwidth=2, binoverlap=10, save_figure=false)
#     # Set bin width and spacing
#     # binwidth = round(nanrange(mdl.Age)/10,sigdigits=1) # Can also set manually, commented out below
#     # binwidth = 100 # Same units as smpl.Age
#     ages = collect(minimum(mdl.Age):binwidth/binoverlap:maximum(mdl.Age))
#     bincenters = ages[1+Int(binoverlap/2):end-Int(binoverlap/2)]
#     spacing = binoverlap

#     # Calculate rates for the stratigraphy of each markov chain step
#     dhdt_dist = Array{Float64}(undef, length(ages)-binoverlap, config.nsteps)
#     @time for i=1:config.nsteps
#         heights = linterp1(reverse(agedist[:,i]), reverse(mdl.Height), ages)
#         dhdt_dist[:,i] .= abs.(heights[1:end-spacing] - heights[spacing+1:end]) ./ binwidth
#     end

#     # Find mean and 1-sigma (68%) CI
#     dhdt = nanmean(dhdt_dist,dim=2)
#     dhdt_50p = nanmedian(dhdt_dist,dim=2)
#     dhdt_025p = nanpctile(dhdt_dist,2.5,dim=2) # Minus 2-sigma (2.5th percentile)
#     dhdt_975p = nanpctile(dhdt_dist,97.5,dim=2) # Plus 2-sigma (97.5th percentile)
#     dhdt_16p = nanpctile(dhdt_dist,15.865,dim=2) # Minus 1-sigma (15.865th percentile)
#     dhdt_84p = nanpctile(dhdt_dist,84.135,dim=2) # Plus 1-sigma (84.135th percentile)
#     # Other confidence intervals (10:10:50)
#     # dhdt_20p = nanpctile(dhdt_dist,20,dim=2)
#     # dhdt_80p = nanpctile(dhdt_dist,80,dim=2)
#     # dhdt_25p = nanpctile(dhdt_dist,25,dim=2)
#     # dhdt_75p = nanpctile(dhdt_dist,75,dim=2)
#     # dhdt_30p = nanpctile(dhdt_dist,30,dim=2)
#     # dhdt_70p = nanpctile(dhdt_dist,70,dim=2)
#     # dhdt_35p = nanpctile(dhdt_dist,35,dim=2)
#     # dhdt_65p = nanpctile(dhdt_dist,65,dim=2)
#     # dhdt_40p = nanpctile(dhdt_dist,40,dim=2)
#     # dhdt_60p = nanpctile(dhdt_dist,60,dim=2)
#     # dhdt_45p = nanpctile(dhdt_dist,45,dim=2)
#     # dhdt_55p = nanpctile(dhdt_dist,55,dim=2)

#     # Plot results
#     hdl = plot(bincenters,dhdt, label="Mean", color=:black, linewidth=2)
#     plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(0,0.5,:darkblue), linealpha=0, label="68% CI")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_025p; reverse(dhdt_975p)], fill=(0,0.2,:darkgrey), linealpha=0, label="95 CI")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_25p; reverse(dhdt_75p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_30p; reverse(dhdt_70p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_35p; reverse(dhdt_65p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_40p; reverse(dhdt_60p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
#     # plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_45p; reverse(dhdt_55p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
#     plot!(hdl,bincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
#     plot!(hdl, xlabel="Age ($(section.Age_Unit))", ylabel="Emplacement Rate", fg_color_legend=:white, framestyle=:box)
#     plot!(hdl, xflip=true)
#     if save_figure
#         savefig(hdl,section.Path*"DepositionRateModel.pdf")
#     end
#     return hdl
# end

function plot_accumulation_model(
    section,
    config,
    mdl,
    agedist;
    binwidth = 2,
    stepfrac = 0.1,   # fraction of binwidth between windows
    save_figure = false
)

    # ------------------------------------------------------------------
    # 1. Define sliding window geometry
    # ------------------------------------------------------------------
    agemin, agemax = extrema(mdl.Age)
    halfw = binwidth / 2
    Δt = binwidth * stepfrac
    println("agemin: $agemin, agemax: $agemax, halfw: $halfw, Δt: $Δt")
    centers = collect(range(agemin + halfw, agemax - halfw, step = Δt))
    println("Number of windows: $(length(centers))")
    if length(centers) < 1
        @warn "Age range too short for accumulation-rate estimation"
        return plot()
    end

    nwin = length(centers)
    nstep = config.nsteps

    # ------------------------------------------------------------------
    # 2. Allocate result array
    # ------------------------------------------------------------------
    dhdt_dist = Array{Float64}(undef, nwin, nstep)

    # ------------------------------------------------------------------
    # 3. Sliding-window rate calculation
    # ------------------------------------------------------------------
    @time for i in 1:nstep
        H = reverse(mdl.Height)
        T = reverse(agedist[:, i])

        for (j, tc) in enumerate(centers)
            t1 = tc - halfw
            t2 = tc + halfw

            h1 = linterp1(T, H, t1)
            h2 = linterp1(T, H, t2)

            dhdt_dist[j, i] = abs(h2 - h1) / binwidth
        end
    end

    # ------------------------------------------------------------------
    # 4. Ensemble statistics
    # ------------------------------------------------------------------
    dhdt_mean = nanmean(dhdt_dist, dims = 2)
    dhdt_med  = nanmedian(dhdt_dist, dims = 2)
    dhdt_16p  = nanpctile(dhdt_dist, 15.865, dims = 2)
    dhdt_84p  = nanpctile(dhdt_dist, 84.135, dims = 2)

    # ------------------------------------------------------------------
    # 5. Plot
    # ------------------------------------------------------------------
    hdl = plot(
        centers,
        dhdt_mean,
        label = "Mean",
        color = :black,
        linewidth = 2
    )

    plot!(
        hdl,
        [centers; reverse(centers)],
        [dhdt_16p; reverse(dhdt_84p)],
        fill = (0, 0.5, :darkblue),
        linealpha = 0,
        label = "68% CI"
    )

    plot!(
        hdl,
        centers,
        dhdt_med,
        label = "Median",
        color = :grey,
        linewidth = 1
    )

    plot!(
        hdl,
        xlabel = "Age ($(section.Age_Unit))",
        ylabel = "Emplacement Rate",
        framestyle = :box,
        fg_color_legend = :white,
        xflip = true
    )

    if save_figure
        savefig(hdl, section.Path * "DepositionRateModel.pdf")
    end

    return hdl
end

function plot_posterior_paths(section, agedist, mdl; n_paths=10, color="grey", linealpha=0.4, save_figure=false)
    # draw n_paths random paths from the posterior distribution with no replacement
    paths = sample(1:size(agedist,2),10,replace=false)
    hdl = plot( agedist[:,paths], mdl.Height, linecolor=color, linealpha=linealpha, label="")
    plot!(hdl, xflip=true)
    plot!(hdl, xlabel="Age ($(section.Age_Unit))", ylabel="Height ($(section.Height_Unit))", fg_color_legend=:white, framestyle=:box)
    if save_figure
        savefig(hdl,section.Path*"PosteriorPaths.pdf")
    end
    return hdl
end
