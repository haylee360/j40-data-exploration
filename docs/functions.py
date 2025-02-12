import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns         
import numpy as np
import os
import libpysal as lps 
from libpysal.weights import W 
from esda.getisord import G_Local
from functools import reduce


def state_gstar(state):
    '''
    state: chr string of state name that you wish to calculate a g star statistic for. Using total burdens and total criteria
    '''
    # Import data
    base_dir = "/capstone/justice40"
    # base_dir = "~/MEDS/justice40/data-exploration"

    # 2.0 communities files (from current CEJST website)
    comm_v2 = pd.read_csv(os.path.join(base_dir, "data", "2.0-communities.csv"))

    # Version 2.0 shapefile data
    v2_geo = gpd.read_file(os.path.join(base_dir, "data", "2.0-shapefile-codebook", "usa", "usa.shp"))

    # Filter to state
    state_df = comm_v2[comm_v2['State/Territory'] == state]

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
    state_geo = v2_geo[v2_geo['SF'] == state]
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


def map_burden(data):
    # Import data
    base_dir = "/capstone/justice40"
    # base_dir = "~/MEDS/justice40/data-exploration"

    # 2.0 communities files (from current CEJST website)
    comm_v2 = pd.read_csv(os.path.join(base_dir, "data", "2.0-communities.csv"))

    # Version 2.0 shapefile data
    v2_geo = gpd.read_file(os.path.join(base_dir, "data", "2.0-shapefile-codebook", "usa", "usa.shp"))

    # Create vector of 48 state names
    state_names = ["Alabama", "Arkansas", "Arizona", "California", "Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Iowa", "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", "Louisiana", "Massachusetts", "Maryland", "Maine", "Michigan", "Minnesota", "Missouri", "Mississippi", "Montana", "North Carolina", "North Dakota", "Nebraska", "New Hampshire", "New Jersey", "New Mexico", "Nevada", "New York", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", "West Virginia", "Wyoming"]

    # Filter to 48 continental states
    comm_states = comm_v2[comm_v2['State/Territory'].isin(state_names)]

    # Climate Change
    df = comm_states[['Expected agricultural loss rate (Natural Hazards Risk Index) (percentile)', 
                'Expected building loss rate (Natural Hazards Risk Index) (percentile)', 
                'Expected population loss rate (Natural Hazards Risk Index) (percentile)', 
                'Share of properties at risk of flood in 30 years (percentile)', 
                'Share of properties at risk of fire in 30 years (percentile)',
                'State/Territory',
                'Census tract 2010 ID']]
    df['mean'] = df[['ag_loss', 'building_loss', 'population_loss', 'flood_risk', 'fire_risk']].mean(axis=1)

    # geometry, GEOID10
    geo_join = v2_geo[['GEOID10', 'geometry']].rename(columns={'GEOID10':'tract_id'})
    geo_join['tract_id'] =geo_join['tract_id'].astype('int64')

    # Merge dfs and reapply geodataframe
    df = pd.merge(df, geo_join, how='left', on='tract_id')
    df = gpd.GeoDataFrame(df)

    # Drop geometry NAs
    df_clean = df[(df.geometry.type == 'Polygon') |( df.geometry.type == 'MultiPolygon')]

    # Remove NAs for G star calc
    df_clean = df_clean.dropna(subset=['mean'])

    # Create weights using Queen method 
    w = lps.weights.Queen(df_clean['geometry'])

    # Run the Getis Ord test for burdens
    g_local = G_Local(df_clean['mean'], w=w, transform='R', permutations=9999)

    # Add statistics to to a new df
    results = df_clean.copy()
    results['standardized_gstar'] = g_local.Zs 
    results['p_norm'] = g_local.p_norm 
    results['p_sim'] = g_local.p_sim

    # PLOT 

    # Define 5 distinct colors
    colors = ["#0000bf", "#4b4bfd", "#efeada", "#e85858", "#a60202"]  # Dark Blue, Blue, Tan, Red, Dark Red

    # Create a discrete colormap
    cmap = mcolors.ListedColormap(colors)

    # Define color boundaries
    # bounds = [0, 1, 2, 3, 4, 5]  # 5 bins
    bounds = [-3, -2.58, -1.96, 1.96, 2.58, 3]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)


    # Initialize figure
    fig, ax = plt.subplots(figsize=(9,5))

    # Remove axis for a cleaner map and set title
    ax.axis('off')
    ax.set_title('Climate Change G Star Heatmap',
                fontsize=14)

    # Plot NY state and color by number of spills 
    plot = results.plot(ax=ax,
                    column='standardized_gstar',
                    cmap=cmap,
                    norm=norm,
                    legend=False)
                    # legend_kwds={
                    #     'label':'Number of Burdens'
                    # })


    # Get the colorbar from the plot
    cbar = plot.get_figure().colorbar(plot.collections[0], ax=ax)

    # # Modify the colorbar tick labels
    cbar.set_ticks([-2.78, -2.27, 0, 2.27, 2.78])
    cbar.set_ticklabels(['Cold Spot 99%+', 'Cold Spot 95%', 'Not Significant', 'Hot Spot 95%', 'Hot Spot 99%+']) 

    plt.show()
