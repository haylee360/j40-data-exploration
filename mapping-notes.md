# Mapping Notes

1. CEJST uses Tippecanoe to generate tiles
    - This tool only works with vector data. 
    - I'm gonna see if there's an easy way to modify their existing tile generation script to just add another layer. 
    - I think this might be more fruitful than trying to mess with the scoring process. We don't actually want to change anything about how the data that goes into the score and how it's calculated, we just want to add functionality on top of what already exists. 

2. Stripped down data frame is still 1 GB. 
    - I stripped down the national data frame to only the geometries and standardized g star for burdens, and the file is still 1 GB as a GeoJSON. Not sure there's really a good way to cut down on size.

3. Didn't get anywhere with making an interactive map elsewhere...
    - Folium ran for 30+ minutes and never got anywhere. Worth trying again, I didn't put that much troubleshooting into it
    - Also worth checking out plotly 

4. There actually are 30 burdens according to the CEJST website
    - I guess when we filter down by state, we only really ever see approx. ~20

5. Need to filter down to just continental U.S. for G star tuning