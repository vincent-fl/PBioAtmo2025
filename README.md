# PBioAtmo2025
Practical master class to learn the practice of eddy covariance

The objectives were to:
* Apply the eddy covariance method in a practical context
* Understand the full processing chain: 
	10 Hz raw data → 30 min fluxes → quality control and u* filtering → gap filled dataset
* Explore relationships between gas fluxes and meteorological variables
* Evaluate suitable gap filling methods for data gaps of different lengths

Measurement details:
* Ecosystem: Grassland, surrounded by agricultural fields, houses & streets
* Footprint: ⇾ Radiation only covers the grassland
⇾ Flux footprint partially covers the surrounding land use
* Canopy: Mowed from 53 cm to 23 cm on 15 May 2025
* Eddy covariance system change: Lifted from 2.20 m to 2.58 m on 19 May 2025
* Measurement period: 09 May 2025 to 16 June 2025 (during growing season)
* 3D ultrasonic anemometer (Gill R3-50) 
* LI-COR enclosed gas analyzer (LI-7200)


## gapfilling 
The gapfilling folder contains the code for analyzing the performance of different gapfilling models.
The code can test different gap sizes with different models, whereby the proportion of data that is filled with gaps and the number of repetitions per gap size and model can be selected. The random seeds are saved and can be reused for the next run. \
It is also possible to set the gaps only at certain times of day, although this still contains potential errors.

## fluxanalysis
The fluxanalysis folder contains various scripts with which, among other things, a standardized file with filters for quality control flags and gap filling by XGBoost was created.
It also contains code for plotting the fluxes and fitting the light response curves.
