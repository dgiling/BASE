BASE v2.3
=========

**NOTE (April 2018): This code is now implemented as an R package 'BASEmetab', which can be installed from github with devtools. See repo (https://github.com/dgiling/BASEmetab) and user guide (https://github.com/dgiling/BASEmetab/blob/master/vignettes/BASEmetab.pdf)**

Code accompanying Grace et al. (2015) Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation). Limnology and Oceanography: Methods, 13, 103–114
http://onlinelibrary.wiley.com/doi/10.1002/lom3.10011/full

PDF user guide included in 'BASE' folder. Download the whole folder as a zip file to get started. 
Contact: Darren Giling (darren.giling@idiv.de) or Mike Grace (Mike.Grace@monash.edu).


Recent updates:

January 2018 (v2.3)
**nb: This is the last update to this code - it is now implemented as a package**
- Simplified file handling (no need to separate files before running)
- Smoothing of PAR and DO
- Automated updating of unconverged chains
- Added better physical constraints on K prior

March 2017 (v2.2)
- Fix to the prior distribution for tau

October 2016 (v2.1)
- Added output of rates on sub-daily timescale

July 2016 (v2.0)
- Changes to model structure following findings of Song et al (2016) L&O:M doi: 10.1002/lom3.10112
- Implementated parallel computing of chains

