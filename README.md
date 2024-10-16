# Model_EVI

Repo to fit a double-log function to annual EVI data, following Beck et al. (2006). Model fitting code is Code adapted by ASL from package `greenbrown` (Forkel and Wutzler, 2015), which is not maintained for the current version of R. We have also tried the model from Elmore et al. (2012), but fits seem to be worse than Beck et al. (2006) for these data.

## Structure

`Code/Model EVI.Rmd` generates model fits from the data files in `Raw_data`. Model fits are then archived in `Output`. Note that raw EVI data from LUMCON are too large to push to GitHub, and are therefore not included in this repo.

## References

Beck, P. S., Atzberger, C., Høgda, K. A., Johansen, B., & Skidmore, A. K. (2006). Improved monitoring of vegetation dynamics at very high latitudes: A new method using MODIS NDVI. Remote sensing of Environment, 100(3), 321-334.

Elmore, A. J., Guinn, S. M., Minsley, B. J. & Richardson, A.D. (2012): Landscape controls on the timing of spring, autumn, and growing season length in mid-Atlantic forests. Global Change Biology 18, 656-674.

Forkel, M. & Wutzler, T. (2015) greenbrown - land surface phenology and trend analysis. A package for the R software. Version 2.2, 2015-04-15, <http://greenbrown.r-forge.r-project.org/>.
