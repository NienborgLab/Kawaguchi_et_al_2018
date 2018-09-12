# Overview
This repository contains Matlab code to replicate results in Kawaguchi et al., 2018.

**"Differentiating between models of perceptual decision-making using pupil-size inferred confidence"**
Katsuhisa Kawaguchi, Stephane Clery, Paria Pourriahi, Lenka Seillier, Ralf Haefner and Hendrikje Nienborg
Journal of Neuroscience 31 August 2018, 0735-18; DOI: https://doi.org/10.1523/JNEUROSCI.0735-18.2018 

To reproduce all the results, first, you need to download our dataset (~ 7GB) from https://figshare.com/articles/Kawaguchi_et_al_2018/7053842.

Once you download the dataset, do the followings.

1. unzip the 'data.zip'
2. put the 'E.mat' into the unzipped 'data' folder
3. clone or download this repository
4. set the path as follows:

Kawaguchi_et_al_2018<br>
&nbsp; &nbsp; &nbsp; &nbsp; ├ data<br>
&nbsp; &nbsp; &nbsp; &nbsp; | &nbsp; &nbsp; ├ figure_storage<br>
&nbsp; &nbsp; &nbsp; &nbsp; | &nbsp; &nbsp; └ struct_storage<br>
&nbsp; &nbsp; &nbsp; &nbsp; ├ analysis_code<br>
&nbsp; &nbsp; &nbsp; &nbsp; └ simulation_code<br>

Please note that you need Matlab 2016b or later version to run the code and generated figures are raw (before formatting).

# Licence
The code in this repository is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License](https://creativecommons.org/licenses/by-nc-sa/3.0/).

Our dataset associated with the paper is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License](https://creativecommons.org/licenses/by-nc-nd/3.0/).

This means you are free to use the data to reproduce and validate our results. However, the materials can only be used for verification unless written permission for another use from us is obtained.

# Organization of the code
There are two main code: 'PS_all.m' and 'Model_all.m'. They reproduce original figures (before formatting) in data anaylsis and model simulation, respectively. These code utilize code stored in '../analysis_code/' and '../simulation_code/'.

To replicate our main findings in model simulations, run 'Model_all.m'. If you have already downloaded our dataset, run 'Model_all(1)' in the matlab command line. Otherwise run 'Model_all(0)', but it would take potentially some days to complete, up to your computer environment, as it starts all the 'grid search' and model fitting from scratch.

To replicate our main findings in data analysis, run 'PS_all.m'. This requires you to have downloaded our datasets in the ../data/ directory. Once you get our data and store them in the ../data/ directory, run 'PS_all('all')' in the matlab command line. This would take some hours to complete.

# Organization of the dataset
To replicate our results, please download our datasets from ... and store it in the ../data/ directory.

## Data matrix ('trmat', 'psmat'):
There are two types of analysis matrices: 'trmat' and 'psmat'.
Both matrices contain relevant information about the animals' behaviors and
stimuli used in the task ('trmat') or preprocessed pupil size time-course ('psmat').

### Preprocessing
The pupil size time-course was preprocessed using either low-pass or band-pass filters.
The type of the preprocessing is indicated in the file-name. Please note that there is no '_lowpass_' file for 'trmat', as it was not used in the paper.

### Animals' names
The part of the file-name '_kiwi' or '_mango' represents each animals' data: <br>
Animal A = kiwi, Animal B = mango.

### Matrix contents
Inside of the matrix, each row represents trial.

The columns contain the following relevant information:

1, session number <br>
2, inferred confidence by pupil size with median split: low (0) and high (1) <br>
3, signal strength with sign <br>
4, three types of available reward size (0, 1, 2) <br>
5, choice (1: far, 0: near) <br>
6, correct (1) or error (0) <br>
7, mean pupil size (mean over 250ms before the stimulus offset) <br>
8, after error (0: no, 1: yes) <br>
9, binary available reward size (0: small, 1: large) <br>
10, completed trial number so far <br>
11, task difficulty (0: hard, 1: intermediate, 2: easy) <br>

('trmat') 12 - 463, presented stimulus sequence <br>
('psmat') 12 - 643, preprocessed pupil size time-course <br>

## Output struct E from Ralf Haefner's sampling model
You can reproduce the Matlab struct E by following the instruction given in [another repository](https://github.com/katsu1110/sampling_decision).  Namely, after cloning the repos, run 'E = S_Experiment(S_Exp_Para('PK-amplitude10'))'. This takes some hours to complete.

## Directories for saving struct and figure
Some of the analysis store additional structs or figures, which are later used in another code. Therefore, directories ../data/struct_storage/ and ../data/figure_storage need to exist.
