import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import colors
import seaborn as sns         
import numpy as np
import os
import libpysal as lps 
from libpysal.weights import W 
from esda.getisord import G_Local
from functools import reduce


def state_gstar(state, df, geo_df):
    """
    state: chr string of state name that you wish to calculate a g star statistic for. Using total burdens and total criteria
    df: chr string data frame of burdens and criteria
    geo_df: chr string geodataframe of geometries for census tracts
    """

    
    # Filter to state
    state_df = df[df['State/Territory'] == state]

    # Narrow down to desired columns and rename
    # Can cut these down depending on what we want...
    state_df = state_df[['Census tract 2010 ID', 'Total threshold criteria exceeded', 'Total categories exceeded', 'Identified as disadvantaged without considering neighbors', 'Identified as disadvantaged based on neighbors and relaxed low income threshold only', 'Identified as disadvantaged', 'Percentage of tract that is disadvantaged by area', 'Share of neighbors that are identified as disadvantaged', 'Total population']].rename(columns={
        'Census tract 2010 ID': 'tract_id',
        'Total threshold criteria exceeded':'total_criteria', 
        'Total categories exceeded':'total_burdens', 
        'Identified as disadvantaged without considering neighbors':'dac_wo_neighbors', 
        'Identified as disadvantaged based on neighbors and relaxed low income threshold only':'dac_relaxed_low_income', 
        'Identified as disadvantaged':'dac', 
        'Percentage of tract that is disadvantaged by area':'percent_area_dac', 
        'Share of neighbors that are identified as disadvantaged':'share_of_dac_neighbors', 
        'Total population':'total_population'
    })

    # state_df = state_df[state_df['dac'] == 1]

    state_df = state_df[state_df['total_population'] > 100]



    # Print Summary stats for burdens and criteria
    print(f"{state} Burden summary stats:\n{state_df['total_burdens'].describe()}")
    print(f"\n{state} Criteria summary stats:\n{state_df['total_criteria'].describe()}")


    # Filter shapefile to state and select desired columns
    state_geo = geo_df[geo_df['SF'] == state]
    state_geo = state_geo[['GEOID10', 'geometry']]

    # Change census tract id type so we can merge later on
    state_geo = state_geo.rename(columns={'GEOID10':'tract_id'})
    state_geo['tract_id'] = state_geo['tract_id'].astype('int64')
    

    # Merge geometry df and burden df
    state_complete = pd.merge(state_df, state_geo, on='tract_id')

    # At some point, we lost our geodataframe type. Add it back so we can map
    # I think this was because of the merge, I just didn't wanna go back and fix it
    state_complete = gpd.GeoDataFrame(state_complete)


    # LISA
    # Drop geometry NAs
    state_clean = state_complete[(state_complete.geometry.type == 'Polygon') |( state_complete.geometry.type == 'MultiPolygon')]

    # Create weights using Queen method 
    w = lps.weights.Queen(state_clean['geometry'])

    # Run the Getis Ord test for burdens
    gstar_burden = G_Local(state_clean['total_burdens'], w, transform='R', permutations=9999) # Set transform=R because the queen weights are row-standardized weights.

    # Run the Getis Ord test for criteria
    gstar_crit = G_Local(state_clean['total_criteria'], w, transform='R', permutations=9999)

    # Add burden stats to to a new df
    results = state_clean.drop('tract_id', axis=1).copy()
    results['standardized_gstar_burd'] = gstar_burden.Zs 
    results['p_norm_burd'] = gstar_burden.p_norm # p-value assuming normal distribution

    # Add criteria stats
    results['standardized_gstar_crit'] = gstar_crit.Zs 
    results['p_norm_crit'] = gstar_crit.p_norm # p-value assuming normal distribution

    #  Results of number of burdens in each bin
    results_95 = results[(abs(results['standardized_gstar_burd']) >=1.96) & (abs(results['standardized_gstar_burd']) < 2.58)]
    results_99 = results[(abs(results['standardized_gstar_burd']) >=2.58)]

    # Results for criteria
    results_95_crit = results[(abs(results['standardized_gstar_crit']) >=1.96) & (abs(results['standardized_gstar_burd']) < 2.58)]
    results_99_crit = results[(abs(results['standardized_gstar_crit']) >=2.58)]

    # Print results
    print(f"\nNumber of tracts in the 95th-99th percentile for burdens: {len(results_95)}")
    print(f"Number of tracts in the 99th+ percentile for burdens: {len(results_99)}")

    print(f"\nNumber of tracts in the 95th-99th percentile for criteria: {len(results_95_crit)}")
    print(f"Number of tracts in the 99th+ percentile for criteria: {len(results_99_crit)}")

    # Plot all four graphs side by side

    # Initialize 2x2 figure grid
    fig, axs = plt.subplots(2, 2, figsize=(9, 5))

    # First plot (top-left)
    axs[0, 0].axis('off')  # Hide axis for a cleaner map
    axs[0, 0].set_title(f"{state} Burdens", fontsize=14)
    results.plot(ax=axs[0, 0], 
                column='total_burdens', 
                cmap='coolwarm', norm=None, legend=True, 
                legend_kwds={'label': 'Number of Burdens'})

    # Second plot (top-right)
    axs[0, 1].axis('off')  # Hide axis for a cleaner map
    axs[0, 1].set_title(f"{state} Criteria", fontsize=14)
    results.plot(ax=axs[0, 1], 
                column='total_criteria', 
                cmap='coolwarm', norm=None, legend=True, 
                legend_kwds={'label': 'Number of Criteria'})


    # Define custom bins for gstar
    bins = [-3, -2.58, -1.96, 0, 1.96, 2.58, 3]

    # Create a colormap
    cmap = plt.get_cmap('coolwarm')

    # Create a norm to map standardized_g_star values to the bins
    norm = colors.BoundaryNorm(boundaries=bins, ncolors=cmap.N)

    # Third plot (bottom-left)
    axs[1, 0].axis('off')  # Hide axis for a cleaner map
    axs[1, 0].set_title(f"{state} Gstar Burden Hotspots", fontsize=14)
    plot1 = results.plot(ax=axs[1, 0], 
                column='standardized_gstar_burd', 
                cmap=cmap, norm=norm, legend=False)
                # legend_kwds={'label': 'Standardized G star'})

    # Get the colorbar from the plot
    cbar = plot1.get_figure().colorbar(plot1.collections[0], ax=axs[1,0])

    # Modify the colorbar tick labels
    cbar.set_ticklabels([' ', '99th+', '95th', 'low', '95th', '99th+', ' ']) 

    # Fourth plot (bottom-right)
    axs[1, 1].axis('off')  # Hide axis for a cleaner map
    axs[1, 1].set_title(f"{state} Gstar Criteria Hotspots", fontsize=14)
    plot2 = results.plot(ax=axs[1, 1], 
                column='standardized_gstar_crit', 
                cmap=cmap, norm=norm, legend=False) 
                # legend_kwds={
                #     'label': 'Standardized G star'
                # })

    # Get the colorbar from the plot
    cbar = plot2.get_figure().colorbar(plot2.collections[0], ax=axs[1,1])

    # Modify the colorbar tick labels
    cbar.set_ticklabels([' ', '99th+', '95th', 'low', '95th', '99th+', ' ']) 

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the final figure
    plt.show()

