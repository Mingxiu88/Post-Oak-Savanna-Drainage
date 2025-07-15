"POS_drainage_GitHub.m" includes 2 sections:
1. estimate drainage using three methods: unit gradient approach, arithmetic mean of K method, and harmonic mean of K method.
2. plot seasonal neutron probe measurement profile, as well as statistical analysis using 3-way ANOVA.

"alignWithDates.m" is a function used in "POS_drainage_GitHub.m" to align data with valid dates.

Data include:
1. daily matric potential (cm) measured at 4 sites (site1: savanna; site2: woodland mosaic under-canopy; site3: woodland mosaic inter-canopy; site4: woodland mosaic under-canopy)
2. daily rainfall data (cm): daily_rain_CR6.csv
3. neutron probe data: NMM measurements_GusEngeling_matlab.xlsx (conversion equation from neutron count ratio to volumetric water content is included in "POS_drainage_GitHub.m")
4. van Genuchten parameters: SoilTexture.xlsx
