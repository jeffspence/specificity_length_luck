# specificity_length_luck
Scripts used for plotting and analyzing data for [Specificity, length and luck drive gene rankings in association studies](https://www.nature.com/articles/s41586-025-09703-7).



To generate the plots, first run `python generate_data.py` (it will take a few hours), then run the R script `plotting_scripts.R` (should take <15 minutes).  `python generate_data.py` will print out filenames that it is generating as it generate them. Before running `plotting_scripts.R` it will be necessary to change the working directory to the repo's directory (on line 19).   Various processed data are available in the `data/` directory.

The required packages are imported at the tops of each of those two scripts.  The python packages with version numbers are listed in `requirements.txt`.  R packages are listed in `r_requirements.txt`. The packages should be generally easy to download and install, except for `fastDTWF` which must be installed as described [here](https://github.com/jeffspence/fastDTWF).  For these analyses we used python v.3.9.23, and R v.4.4.0. 
