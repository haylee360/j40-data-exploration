def map_burden(burden):
    """
    This function accepts one of the following burden names:
    'Climate Change', 'Energy', 'Health', 'Housing', 'Legacy Pollution', 'Transportation', 'Waste and Wastewater', or 'Workforce Development'.
    
    Arguments:
    burden_name (str): The name of the burden, which must be one of the keys in the burden_names dictionary.
    
    Raises:
    ValueError: If the burden_name is not one of the valid options.
    """
    burden_names = {
        'Climate Change':'cc',
        'Energy':'energy',
        'Health':'health',
        'Housing':'housing',
        'Legacy Pollution':'lp',
        'Transportation':'transport',
        'Waste and Wastewater':'ww',
        'Workforce Development':'wd'
    }

    # Check if the input is a valid key in the dictionary
    if burden not in burden_names:
        raise ValueError(f"Invalid burden name: '{burden}'. Must be one of: {', '.join(burden_names.keys())}.")
    
    # If valid, process the corresponding value
    print(f"Processing {burden} burden with code {burden_names[burden]}")


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

    # Need to pull the column names from the data frame and use them to subset comm_states. 

    

    # Don't need to give a data frame, give one of 8 selections and store the dfs in the function, then select the appropraite df based on the chr string

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
