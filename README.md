## --< EXCALIBUR >--

EXCALIBUR stands for EXoplanet CALIbration Bayesian Unified Retrieval.

The EXCALIBUR pipeline enables comparative exoplanet studies by uniformly processing input catalogs of observations to create science data products. EXCALIBUR preserves the chain of inference with persistent data products for each processing step that are tagged to a unique identifier linking them to a specific compute instance and GitHub change set. EXCALIBUR catalog releases are available on the EXCALIBUR portal hosted at IPAC https://excalibur.ipac.caltech.edu. 

EXCALIBUR is an event driven pipeline where the events are defined as changes in data or algorithms; when events are detected, dependencies affected by the changes are re-processed. This makes EXCALIBUR a useful tool for sensitivity testing at catalog scales. EXCALIBUR includes CERBERUS, a line-by-line, spherical-geometry radiative transfer code modeling exoplanet atmospheres and a Bayesian parameter retrieval and model selection package. EXCALIBUR also includes a UI that provides visualization of data products. 

EXCALIBUR papers include:
- https://ui.adsabs.harvard.edu/abs/2024MNRAS.531...35M/abstract
- https://ui.adsabs.harvard.edu/abs/2022ApJ...941L...5E/abstract
- https://ui.adsabs.harvard.edu/abs/2022AJ....163...22H/abstract
- https://ui.adsabs.harvard.edu/abs/2021AJ....162...91E/abstract
- https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.3316B/abstract
- https://ui.adsabs.harvard.edu/abs/2021AJ....161..213S/abstract

Current EXCALIBUR processing capabilities include transit spectroscopy and phase curves for selected instruments. EXCALIBUR is under active development with ongoing capability enhancements, and we welcome contributors and collaborations. 


### --< become a member >--

1. log into https://github-fn.jpl.nasa.gov (may require a lot of steps to get there)
2. ask to become a member of the esp worker team
3. clone this repository

#### --< security certificate >--

If you're going to run a personal pipeline, you'll probably need a certificate.
You definitely need it if you're going to command the main pipeline to run tasks.

1. Create certificate using a script in the github repository
  In esp directory, run `bash/make_cert.sh ${HOME}/.ssh/my_excalibur_identity.pem`
  This creates a file in your .ssh folder, and also in /proj/sdp/data/certs/
2. Convert your .pem to .p12 file (might not always be necessary, but can't hurt)
  Run 'openssl pkcs12 -export -in my_excalibur_identity.pem -out my_excalibur_identity.p12'
  It will ask for a password that you use later when loading the certificate.
3. Copy the .p12 file to your home machine, using scp.
4. Import the .p12 file into your browser
 For Firefox, go to Settings / Security. Certificates are down near the bottom.
 For Safari, you import it via Keychain Access (in /Applications/Utilities)
 For Chrome, requires .p12 file.

Note that for commanding the main pipeline, you need to use port 8085 instead of port 8080 (in the excalibur URL).

### --< Organization >--

#### --< Data >--

Some of the items we do not maintain (external) are static. Static executables and libraries are kept out of the code management.

- /proj/sdp/bin (external executables such as pymc, ...)
- /proj/sdp/lib (external libraries)
- /proj/sdp/pkg (external source code and build area)
- /proj/sdp/data/cal (instrument reference calibration files such as STSCI files, ...)
- /proj/sdp/data/res (resource files such as web interface tools, ...)
- /proj/sdp/data/sci (on-disk data to be processed such as private datasets, ...)

#### --< TASKS / algorithms >--

RUNTIME
- create (loads/sets up the overall parameter file for excalibur, e.g. which filters to run)
- autofill (applies the overall parameter file to each target)

TARGET
- create (IDs and filters)
- autofill (prior system information from NEXSCI)
- scrape (download and save available data)

SYSTEM
- validate (checks for system parameters completeness)
- finalize (delivers a comprehensive set of system parameters uniformly formatted, allows parameters over-ride)

ANCILLARY
- estimate (calculates a set of derived parameters, e.g. the scale height)

ARIEL
- sim_spectrum (simulates Ariel observations)

DATA
- collect (sort data according to filters)
- calibration (extraction, wavelength solution, noise assessment)
- timing (transit, eclipse, full phase curve detection)

TRANSIT
- normalization (scaling of stellar spectrum to Out Of Transit relative quantities)
- whitelight (orbital solution and instrumental behavior recovery)
- spectrum (exoplanet transmission spectrum recovery)

ECLIPSE
- normalization (scaling of stellar spectrum to Out Of Transit relative quantities)
- whitelight (orbital solution and instrumental behavior recovery)
- spectrum (exoplanet emission spectrum recovery)

PHASECURVE
- normalization (scaling of stellar spectrum to Out Of Transit relative quantities)
- whitelight (orbital solution and instrumental behavior recovery)

CLASSIFIER
-flags (calculate various performance metrics, to help weed out bad data/results)

CERBERUS
- xslib (cross section library from EXOMOL and HITEMP/HITRAN)
- atmos (model selection and atmospheric content recovery)
- results (plots the results of the atmos retrieval)
- analysis (plots mass-metallicity and other summary plots for the whole population)

### --< Source Code >--

[Github](https://github-fn.jpl.nasa.gov/EXCALIBUR/esp)
.
