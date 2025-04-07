# EXCALIBUR
## Welcome
EXCALIBUR stands for EXoplanet CALIbration Bayesian Unified Retrieval.

The EXCALIBUR pipeline enables comparative exoplanet studies by uniformly processing input catalogs of observations to create science data products. EXCALIBUR preserves the chain of inference with persistent data products for each processing step that are tagged to a unique identifier linking them to a specific compute instance and GitHub change set. EXCALIBUR catalog releases are available on the [EXCALIBUR portal](https://excalibur.ipac.caltech.edu) hosted at IPAC. 

EXCALIBUR is an event driven pipeline where the events are defined as changes in data or algorithms; when events are detected, dependencies affected by the changes are re-processed. This makes EXCALIBUR a useful tool for sensitivity testing at catalog scales. EXCALIBUR includes CERBERUS, a line-by-line, spherical-geometry radiative transfer code modeling exoplanet atmospheres and a Bayesian parameter retrieval and model selection package. EXCALIBUR also includes a UI that provides visualization of data products. 

EXCALIBUR papers include:
- [Comparing transit spectroscopy pipelines at the catalogue level: evidence for systematic differences (Mugnai et al 2024, MNRAS)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.531...35M/abstract)
- [A Temperature Trend for Clouds and Hazes in Exoplanet Atmospheres (Estrela et al 2022, ApJL)](https://ui.adsabs.harvard.edu/abs/2022ApJ...941L...5E/abstract)
- [Characterization of an Instrument Model for Exoplanet Transit Spectrum Estimation through Wide-scale Analysis on HST Data (Huber-Feely et al. 2022)](https://ui.adsabs.harvard.edu/abs/2022AJ....163...22H/abstract)
- [Detection of Aerosols at Microbar Pressures in an Exoplanet Atmosphere, Estrela et al 2021, AJ](https://ui.adsabs.harvard.edu/abs/2021AJ....162...91E/abstract)
- [Disequilibrium Chemistry in Exoplanet Atmospheres Observed with the Hubble Space Telescope (Roudier et al. 2021 AJ)](https://ui.adsabs.harvard.edu/abs/2021AJ....162...37R/abstract)
- [Detection of an Atmosphere on a Rocky Exoplanet (Swain et al 2021, AJ)](https://ui.adsabs.harvard.edu/abs/2021AJ....161..213S/abstract)

Current EXCALIBUR processing capabilities include transit spectroscopy and phase curves for selected instruments. EXCALIBUR is under active development with ongoing capability enhancements, and we welcome contributors and collaborations.

## Table of Content
[Private Pipeline Disk Mapping](disk_mapping)
