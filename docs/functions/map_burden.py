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

def map_burden(burden, df, geo_df):
    """
    This function accepts one of the following burden names:
    'Climate Change', 'Energy', 'Health', 'Housing', 'Legacy Pollution', 'Transportation', 'Waste and Wastewater', or 'Workforce Development'.
    
    Arguments:
    burden (str): The name of the burden, which must be one of the keys in the burden_names dictionary.
    df: data frame of burden and criteria
    geo_df: geo data frame of census tract geometry
    
    Raises:
    ValueError: If the burden_name is not one of the valid options.
    """


    # Create vector of 48 state names
    state_names = ["Alabama", "Arkansas", "Arizona", "California", "Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Iowa", "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", "Louisiana", "Massachusetts", "Maryland", "Maine", "Michigan", "Minnesota", "Missouri", "Mississippi", "Montana", "North Carolina", "North Dakota", "Nebraska", "New Hampshire", "New Jersey", "New Mexico", "Nevada", "New York", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", "West Virginia", "Wyoming"]

    # Filter to 48 continental states
    comm_states = df[df['State/Territory'].isin(state_names)]

    # Need to pull the column names from the data frame and use them to subset comm_states. 

    

    # Don't need to give a data frame, give one of 8 selections and store the dfs in the function, then select the appropraite df based on the chr string

    # Climate Change
    cc = comm_states[['Expected agricultural loss rate (Natural Hazards Risk Index) (percentile)', 
                'Expected building loss rate (Natural Hazards Risk Index) (percentile)', 
                'Expected population loss rate (Natural Hazards Risk Index) (percentile)', 
                'Share of properties at risk of flood in 30 years (percentile)', 
                'Share of properties at risk of fire in 30 years (percentile)',
                'State/Territory',
                'Census tract 2010 ID']]
    cc =cc.rename(columns={
        'Expected agricultural loss rate (Natural Hazards Risk Index) (percentile)':'ag_loss', 
        'Expected building loss rate (Natural Hazards Risk Index) (percentile)':'building_loss', 
        'Expected population loss rate (Natural Hazards Risk Index) (percentile)':'population_loss', 
        'Share of properties at risk of flood in 30 years (percentile)':'flood_risk', 
        'Share of properties at risk of fire in 30 years (percentile)':'fire_risk',
        'Census tract 2010 ID':'tract_id',
        'State/Territory':'state'
    })
    cc['mean'] = cc[['ag_loss', 'building_loss', 'population_loss', 'flood_risk', 'fire_risk']].mean(axis=1)

    # Energy
    energy = comm_states[['Energy burden (percentile)', 
                        'PM2.5 in the air (percentile)', 
                        'Census tract 2010 ID']]
    energy = energy.rename(columns={
        'Energy burden (percentile)':'energy_burden', 
        'PM2.5 in the air (percentile)':'pm_25',
        'Census tract 2010 ID':'tract_id'
        })
    energy['mean'] = energy[['energy_burden', 'pm_25']].mean(axis=1)

    # Housing
    housing = comm_states[['Housing burden (percent) (percentile)', 
                    'Share of homes with no kitchen or indoor plumbing (percentile)', 
                    'Percent pre-1960s housing (lead paint indicator) (percentile)',
                    'Census tract 2010 ID']]
    housing = housing.rename(columns={
        'Housing burden (percent) (percentile)':'housing_burden',
        'Share of homes with no kitchen or indoor plumbing (percentile)':'no_plumbing', 
        'Percent pre-1960s housing (lead paint indicator) (percentile)':'lead_paint',
        'Census tract 2010 ID':'tract_id'
    })
    housing['mean'] = housing[['housing_burden', 'no_plumbing', 'lead_paint']].mean(axis=1)

    # Health
    health = comm_states[['Current asthma among adults aged greater than or equal to 18 years (percentile)', 
                    'Diagnosed diabetes among adults aged greater than or equal to 18 years (percentile)', 
                    'Coronary heart disease among adults aged greater than or equal to 18 years (percentile)',
                    'Low life expectancy (percentile)',
                    'Census tract 2010 ID']]
    health = health.rename(columns={
        'Current asthma among adults aged greater than or equal to 18 years (percentile)':'asthma', 
        'Diagnosed diabetes among adults aged greater than or equal to 18 years (percentile)':'diabetes', 
        'Coronary heart disease among adults aged greater than or equal to 18 years (percentile)':'heart_disease',
        'Low life expectancy (percentile)':'low_life_expectancy',
        'Census tract 2010 ID':'tract_id'
    })
    health['mean'] = health[['asthma', 'diabetes', 'heart_disease', 'low_life_expectancy']].mean(axis=1)

    # Legacy Pollution
    lp = comm_states[['Is there at least one abandoned mine in this census tract?', 
                'Is there at least one Formerly Used Defense Site (FUDS) in the tract?', 
                'Proximity to hazardous waste sites (percentile)', 
                'Proximity to NPL (Superfund) sites (percentile)', 
                'Proximity to Risk Management Plan (RMP) facilities (percentile)',
                'Census tract 2010 ID']]
    lp = lp.rename(columns={
        'Is there at least one abandoned mine in this census tract?':'abandoned_mines', 
        'Is there at least one Formerly Used Defense Site (FUDS) in the tract?':'defense_site', 
        'Proximity to hazardous waste sites (percentile)':'hazardous_waste', 
        'Proximity to NPL (Superfund) sites (percentile)':'superfund_sites', 
        'Proximity to Risk Management Plan (RMP) facilities (percentile)':'rmp_facilites',
        'Census tract 2010 ID':'tract_id'
    })
    lp['mean'] = lp[['hazardous_waste', 'superfund_sites', 'rmp_facilites']].mean(axis=1)

    # Transportation
    transport = comm_states[['Diesel particulate matter exposure (percentile)', 
                        'DOT Travel Barriers Score (percentile)', 
                        'Traffic proximity and volume (percentile)',
                        'Census tract 2010 ID']]
    transport = transport.rename(columns={
        'Diesel particulate matter exposure (percentile)':'diesel_pm', 
        'DOT Travel Barriers Score (percentile)':'travel_barriers', 
        'Traffic proximity and volume (percentile)':'traffic_proximity',
        'Census tract 2010 ID':'tract_id'
    })
    transport['mean'] = transport[['diesel_pm', 'travel_barriers', 'traffic_proximity']].mean(axis=1)

    # Waste and Wastewater
    ww = comm_states[['Leaky underground storage tanks (percentile)', 'Wastewater discharge (percentile)', 'Census tract 2010 ID']]
    ww =ww.rename(columns={
        'Leaky underground storage tanks (percentile)':'leaky_storage_tanks', 
        'Wastewater discharge (percentile)':'wastewater_discharge',
        'Census tract 2010 ID':'tract_id'
    })
    ww['mean'] = ww[['leaky_storage_tanks', 'wastewater_discharge']].mean(axis=1)

    # Workforce Development
    wd = comm_states[['Linguistic isolation (percent) (percentile)', 
                'Low median household income as a percent of area median income (percentile)', 
                'Percent of individuals below 200% Federal Poverty Line (percentile)', 
                'Unemployment (percent) (percentile)',
                'Census tract 2010 ID']]
    wd = wd.rename(columns={
        'Linguistic isolation (percent) (percentile)':'ling_isolation', 
        'Low median household income as a percent of area median income (percentile)':'low_income', 
        'Percent of individuals below 200% Federal Poverty Line (percentile)':'poverty', 
        'Unemployment (percent) (percentile)':'unemployment',
        'Census tract 2010 ID':'tract_id'
    })
    wd['mean'] = wd[['ling_isolation', 'low_income', 'poverty', 'unemployment']].mean(axis=1)

    burden_names = {
        'Climate Change': cc,
        'Energy': energy,
        'Health': health,
        'Housing': housing,
        'Legacy Pollution': lp,
        'Transportation': transport,
        'Waste and Wastewater': ww,
        'Workforce Development': wd
    }

    # Check if the input is a valid key in the dictionary
    if burden not in burden_names:
        raise ValueError(f"Invalid burden name: '{burden}'. Must be one of: {', '.join(burden_names.keys())}.")
    
    # If valid, process the corresponding value
    print(f"Processing {burden} burden...")

    # Get the appropriate data frame
    data = burden_names[burden]

    # geometry, GEOID10
    geo_join = geo_df[['GEOID10', 'geometry']].rename(columns={'GEOID10':'tract_id'})
    geo_join['tract_id'] =geo_join['tract_id'].astype('int64')

    # Merge dfs and reapply geodataframe
    data = pd.merge(data, geo_join, how='left', on='tract_id')
    data = gpd.GeoDataFrame(data)

    # Drop geometry NAs
    data_clean = data[(data.geometry.type == 'Polygon') |( data.geometry.type == 'MultiPolygon')]

    # Remove NAs for G star calc
    data_clean = data_clean.dropna(subset=['mean'])

    # Create weights using Queen method 
    w = lps.weights.Queen(data_clean['geometry'])

    # Run the Getis Ord test for burdens
    g_local = G_Local(data_clean['mean'], w=w, transform='R', permutations=9999)

    # Add statistics to to a new df
    results = data_clean.copy()
    results['standardized_gstar'] = g_local.Zs 
    results['p_norm'] = g_local.p_norm 
    results['p_sim'] = g_local.p_sim

    # PLOT 

    # Define 5 distinct colors
    colors = ["#0000bf", "#4b4bfd", "#efeada", "#e85858", "#a60202"]  # Dark Blue, Blue, Tan, Red, Dark Red

    # Create a discrete colormap
    cmap = mcolors.ListedColormap(colors)

    # Define color boundaries
    bounds = [-3, -2.58, -1.96, 1.96, 2.58, 3]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Initialize figure
    fig, ax = plt.subplots(figsize=(9,5))

    # Remove axis for a cleaner map and set title
    ax.axis('off')
    ax.set_title(f"{burden} G Star Heatmap",
                fontsize=14)

    # Plot NY state and color by number of spills 
    plot = results.plot(ax=ax,
                    column='standardized_gstar',
                    cmap=cmap,
                    norm=norm,
                    legend=False)

    # Get the colorbar from the plot
    cbar = plot.get_figure().colorbar(plot.collections[0], ax=ax)

    # # Modify the colorbar tick labels
    cbar.set_ticks([-2.78, -2.27, 0, 2.27, 2.78])
    cbar.set_ticklabels(['Cold Spot 99%+', 'Cold Spot 95%', 'Not Significant', 'Hot Spot 95%', 'Hot Spot 99%+']) 

    plt.show()
