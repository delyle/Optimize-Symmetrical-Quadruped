# OptSymQuad

Dependencies:

* [GPOPS-II](http://www.gpops2.com/)
* [SNOPT](https://ampl.com/products/solvers/solvers-we-sell/snopt/)

[Preprint](https://www.biorxiv.org/content/10.1101/2020.04.24.060509v1)

This software uses through-contact optimization (the general approach is described in [this PLOS Comp Biol paper](https://dx.doi.org/10.1371%2Fjournal.pcbi.1007444)) to find symmetrical gaits that optimize a work-based objective. This particular code was used to explore the effect of speed and moment of inertia on optimal gait, as described in [this preprint](https://www.biorxiv.org/content/10.1101/2020.04.24.060509v1). Optimal solutions and their beat number for a large set of moment of inertia and speed are also provided as MATLAB binaries.

This code runs in MATLAB and has been tested on MATLAB 2019a for macOS and MATLAB 2020b for Windows.

## Running the code

1. In MATLAB, `cd` to the OptSymQuad directory
2. `run add_necessary_paths`. This will add the relevant functions to your MATLAB path.
3. You can copy the contents of `add_necessary_paths` to your `startup.m` file to add the paths on MATLAB startup
4. Make sure that GPOPS-II and SNOPT are in your path as well

### Perform your own optimization


#### Run optimization
To run an optimization yourself, first run `MAIN_SolveSymQuad`. Parameters can be adjusted in the first section of the script. With default parameters, the code generates a sparse dataset in a reasonable amount of time (though it may take up to an hour depending on hardware).

#### Compile results
To compile results, run `CompileOptimalSolutions`. This makes a single MATLAB binary with pseudoglobal optimal solutions for easy access. With default parameters in `MAIN_SolveSymQuad`, the script will run correctly as-is. Otherwise, be sure it is pointing to the correct directory (see user parameters in the script)

#### Detect Gait Types
Next, run `MAIN_DetectGaitTypes`. As written, the script uses the gait detection parameters from the paper, and points to a compiled datafile with data from the paper.

This file will map an optimal solution to a gait, based on its beat number and whether it is a "walk" or "run".

If you have generated your own data, be sure to set `data_path` to the correct data file, which should be something like `Data/mf*/BestGaitSolutions*.mat`, where `*` can be different numbers depending on parameters and time that files were generated. Note that the parent data directory here is `Data/mf*`.

If the script runs correctly, it will save a file in the parent data directory called `GaitTypeData*.mat`.

#### Gait Type figure
Finally, run `MAIN_GaitTypeFigure`. As written, the file retrieves data from the paper. If using your own data, make sure `data_path` points to the `GaitTypeData*.mat` file you generated above. Adjust `data_path` in the user inputs section of the script.

This file will generate a figure similar to Fig 3a in the paper.

## MAIN_SolveSymQuad

This script iterates through some combination of normalized Murphy number ($I/mL^2$) and T' ($T\sqrt{g/(2L)}$) to find the optimal gait for each combination.

Default values for the model used in the paper are included
in the USER INPUTS section. Out of the box, T and I are specified over a sparse grid, and n (the number of guesses per T and I combination) is set to a low value. The code can run in a reasonable amount of time. For producing the data in the paper, the code was run over a much denser grid, and 50-150 guesses were used (depending on I and T).

### Running in Parallel
The script can be run with multiple instances of matlab in parallel. It will automatically detect whether the current guess has been started and, if so, will move on to the next guess.

Because of this feature, note that if the script is terminated early andrerun, it will start at the next guess after the one in which it was terminated (even if GPOPS-II did not complete its operation for that guess).

If the script is terminated early, a temporary directory `tmp*`with a (potentially very large) SNOPT text file will remain in pdir. This directory can be safely deleted manually if the script is not running.

## MAIN_DetectGaitTypes

This script uses data collected by CompileOptimalSolutions and determines the gait used for each case.

Data used in the paper are included out of the box as 
`Data/BestGaitSolutions.mat`

If the user wants to compile their own data, run `CompileOptimalSolutions` and then specify the path to the resulting data file on line 15 of `MAIN_DetectGaitTypes`

## MAIN_GaitTypeFigure

This script creates a plot showing gait regions as a function of Murphy number and normalized speed.

Out of the box, it should reproduce figure 3a in the paper. A path to custom data, produced by `MAIN_DetectGaitTypes`, can be specified on line 11 in `MAIN_GaitTypeFigure`.


