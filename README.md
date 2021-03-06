# CS_inference
## To use our scripts, please read README first.

### 1. CS_inference provides all the Matlab scripts and some example files that use compressed sensing algorithms to infer the original singe cell expression matrix in this article.

### 2. /regression_model contains two core scripts for inferring original data by the two regression models.<br>
Among them, <br>
* `l1eq_pd.m` is the Basis Pursuit model, quoted from:<br>
>Candes, EJ, and J. Romberg. "Toolbox ℓ1-MAGIC." California Inst. of Technol., Pasadena, CA (http://www.acm.caltech.edu/l1magic/);<br>
* `Ridge_direct.m` is the Ridge Regression model, quoted from:<br>
>W. Pan, Y. Yuan, J. Goncalves, and G.-B. Stan, A Sparse Bayesian Approach to the Iden- tification of Nonlinear State-Space Systems, IEEE Transaction on Automatic Control, 2015 (to appear). arXiv:1408.3549 (http://arxiv.org/abs/1408.3549) <br>
* These two core calculation scripts are used in both 3 and 4 process below.

### 3. /computer_simulation contains scripts used in computer simulation experiments.<br>
Among them, <br>
* `l1eq_ridge_comparison.m` is the main function script to be executed, and its calculation process is:<br>
>* Step1:  Randomly generate a measurement matrix, and get the corresponding observation matrix;<br>
>* Step2:  Use the two regression models to execute the compressed sensing algorithm to obtain the inferred recover matrix (some other parameters can be set here); <br>
>* Step3:  Calculate the Pearson correlation between the obtained recover matrix and the original singe cell expression matrix, and calculate the average/median/SIR value, and put the results of the two methods in one data matrix.<br>
>* (The implementation of step2 requires `l1eq_simulation_full.m` and `ridge_simulation_full.m` scripts.)<br>
>* See the annotations in the script for input and output details.<br>
* For the inference of a large-scale data set, use `parallel.m` to divide the original matrix into several sub-blocks first and then conduct the restoration (need `ridge_simulation_full.m` script).

### 4. /experiment_inference contains the scripts and a sample file used in the in vitro experiment.<br>
Among them, <br>
* `phi_normalization.mat` is an example file of the normalized measurement matrix (the method of generating this matrix can be found in the "Read alignments and gene-expression estimation" section of the article);<br>
* `l1eq_experiment.m` is the main function script of the Basis Pursuit model to be executed;<br>
* `ridge_experiment.m` is the main function script of the Ridge Regression model to be executed; <br>
>The two functions’ calculation process is: through the provided pool expression matrix and the normalized measurement matrix, use the corresponding regression model to execute the compressed sensing algorithm and infer the original single cell expression matrix.<br>
>* See the annotations in the script for input and output details.


### 5. To make a basic run, please execute a simple `run.m` as (just for an example): <br>
```Matlab
ridge_experiment(input, inputphi, outdir, core)
