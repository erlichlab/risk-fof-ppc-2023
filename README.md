This is the repository for code and data to generate results and analysis for the following paper:

"The rat frontal orienting field dynamically encodes value for economic decisions under risk". 
Bao, C.†, Zhu, X.†, Moller-Mara, J., Li, J., Dubroqua, S., & Erlich, J. C. (2023). https://doi.org/10.1038/s41593-023-01461-x

† These authors contributed equally to the work. 



## Repository contents
The repository comprises: 

1. Behavioral data (including muscimol and optogenetic perturbation data) in the CSV files (in `CSV` folder).
2. Electrophysiological data in Matlab ".mat" files (in `matlab` folder). 
3. R code for generating the plots and statistics for the behavior and the three-agent model (in `R` folder)
4. Matlab (R2022a) code and analysis for electrophysiological recording (in `matlab` folder). 
5. Julia code to generate the plots for the dynamical model, alternative dynamical models and pseudopopulation decoding of lottery magnitude from FOF electrophysiological data (in `julia` folder).

## Instructions for replicating plots and results

Clone this reposistory
```
cd 
git clone https://github.com/erlichlab/risk-fof-ppc-2023
```

### Generating the behavior and the three-agent model related analysis and plots

Download the model fits from the following link: https://figshare.com/s/bce33318aac210a8279c

Move the `brms_model_fits.RData` file to the `risk-fof-ppc-2023` root directory

Create the R enviroment using conda:

First install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba). (Note: mamba isn't required but it is faster than using the base conda setup.)

then in a terminal run:
```
cd  ~/risk-fof-ppc-2023
mamba install -f environment.yml
```

In your R enviroment, set your working directory: `setwd('~/risk-fof-ppc-2023')`.

Load all the library and source code for generating the plots and analysis: `source('R/init.R')`

Then, generate all the R plots for the paper: `source('R/main.R')`

## Generate the electrophysiology analysis and plots

Install [elutils](https://github.com/erlichlab/elutils) and add to the path
```
bash> git clone  https://github.com/erlichlab/elutils.git
matlab> addpath ~/elutils
matlab> cd risk-fof-pcc-2023/matlab
matlab> main_ephys_plots
matlab> lottery_sound_control
```

## Generate the analysis for pseudopopulation decoding of lottery magnitude from FOF ephys data

If you don't have julia, then first install julia (recommended way is to use [juliaup](https://github.com/JuliaLang/juliaup)).

Go into the `julia` folder of this project and run julia like this:

`julia --project=. -tauto` 

```
import Pkg
Pkg.instantiate() # This gets the dependcies
using Pluto
Pluto.run(notebook="fof_lottery_decode.jl")
```

Note: each Pluto notebook is its own environment. The first time you run the notebook, it will download and precompile a lot of dependecies, so it can take quite a while.

## Generate the dynamic model related plots

Follow the instructions above but run  `Pluto.run(notebook="julia/dynamical.jl")`

