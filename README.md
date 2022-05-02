# CRISPRtracrRNA: Robust approach for tracrRNA detection
The CRISPR-Cas9 system is a type II CRISPR system that has rapidly become the most versatile and widespread tool for genome engineering. It consists of two components, the Cas9 effector protein, and a single guide RNA that combines the spacer (for identifying the target) with the tracrRNA, a trans-activating small RNA required for both crRNA maturation and interference. While there are well-established methods for screening Cas effector proteins and CRISPR arrays, the detection of tracrRNA remains the bottleneck in detecting Class 2 systems.
Results: We introduce a new pipeline CRISPRtracer for screening and evaluation of tracrRNA candidates in genomes. This pipeline combines evidence from different components of the Cas9-sgRNA complex. The core is a newly developed structural model via covariance models from sequence-structure alignment of experimentally validated tracrRNAs. As additional evidence, we determine the terminator signal (required for the tracrRNA transcription) and the RNA-RNA interaction between the CRISPR array repeat and the 5’-tail of the tracrRNA. Repeats are detected via an ML-based approach (CRISPRidenifier). As additional evidence, we detect the cassette containing the Cas9 (type II CRISPR systems) and Cas12 (type V CRISPR systems) effector protein. Our tool is the first for detecting tracrRNA for type V systems.


## Getting Started with CRISPRtracrRNA

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

First you need to install Miniconda
Then create an environment and install the required libraries in it


### Creating a Miniconda environment

First you need to install Miniconda for python 3.
Miniconda can be downloaded from here:

https://docs.conda.io/en/latest/miniconda.html

Then Miniconda should be installed. On a linux machine the command is similar to this one:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Then you have to create the CRISPRtracrRNA. The necessary setup is provided in the "environment.yml" file inside the "for_environment" directory

In order to install the corresponding environment one can execute the following command from the "for_environment" directory

```
conda env create -f environment.yml
```

### Additional preparations

CRISPRtracrRNA utilizes CRISPRidentify for CRISPR-array search and CRISPRcasIdentifier for the detection of the cas genes.

You can find the CRISPRidentify tool and its description [here](https://github.com/BackofenLab/CRISPRidentify)

Please make sure that after you downloaded CRISPRidentify its relative path is:

```
tools/CRISPRidentify/CRISPRidentify/CRISPRidentify.py
```



You can find the CRISPRcasIdentifier tool and its description [here](https://github.com/BackofenLab/CRISPRcasIdentifier)

You need to make two steps:

Firstly, you need to download the CRISPRcasIdentifier tool:
```
wget https://github.com/BackofenLab/CRISPRcasIdentifier/archive/v1.1.0.tar.gz
tar -xzf v1.1.0.tar.gz
```

Please make sure that after you downloaded CRISPRcasIdentifier its relative path is:

```
tools/CRISPRcasIdentifier/CRISPRcasIdentifier/CRISPRcasIdentifier.py
```
Secondly, you need to download the models:

Due to GitHub's file size constraints, authors made their HMM and ML models available in Google Drive. You can download them [here](https://drive.google.com/file/d/1YbTxkn9KuJP2D7U1-6kL1Yimu_4RqSl1/view?usp=sharing) and [here](https://drive.google.com/file/d/1Nc5o6QVB6QxMxpQjmLQcbwQwkRLk-thM/view?usp=sharing). Save both tar.gz files inside CRISPRcasIdentifier's directory.

### Activation of the environment

Before running CRISPRtracrRNA you need to activate the corresponding environment.

```
conda activate crispr_tracr_rna_env
```

## Running CRISPRtracrRNA

We prepared the test folder which can be used for the test run.

Example of running CRISPRtracrRNA over a folder of files:

```
python CRISPRtracrRNA.py --input_folder test_input_complete
```

### Flags

You can see the help by using the `-h` option

```

python CRISPRtracrRNA.py -h

```


#### Mandatory flags
The only mandatory parameter which has to be specified is the input folder.

##### Input as a folder of fasta files

* `--input_folder <path_to_the_folder>`


#### Optional flags

##### Output

* `--output_folder [paht_to_the_result_folder]`

Specifies the path and name of the folder with the output results. If not specified the results will appear in "CRISPRtracrRNA_output" folder

For example:
Specifies the path and name of the folder with the output results.
```
python CRISPRtracrRNA.py --input_folder test_input_complete --output_folder CRISPRtracrRNA_output
```

##### Temp folder

* `--temp_folder_path [paht_to_the_result_folder]`

Specifies the path and name of the folder with the intermediate results. If not specified the results will appear in "temp" folder

```
python CRISPRtracrRNA.py --input_folder test_input_complete --temp_folder_path temp
```

##### Run type

* `--run_type [execution_type]`

Specifies the way CRISPRtracrRNA will be executed. The following options are available: "complete_run" and "model_run"

```
python CRISPRtracrRNA.py --input_folder test_input_complete --run_type complete_run
```

```
python CRISPRtracrRNA.py --input_folder test_input_model_only --run_type model_run
```

##### Model type

* `--model_type [model]`

Specifies the type of pretrained model. The following options are available: "II" and "V"

```
python CRISPRtracrRNA.py --input_folder test_input_complete --model_type II
```

##### Anti-repeat similarity threshold

* `--anti_repeat_similarity_threshold [value]`

Specifies the threshold for the similarity between anti-repeat and repeat. The default value is 0.7.

```
python CRISPRtracrRNA.py --input_folder test_input_complete --anti_repeat_similarity_threshold 0.7
```

##### Anti-repeat coverage threshold

* `--anti_repeat_coverage_threshold [value]`

Specifies the threshold for the coverage of anti-repeat. The default value is 0.6.

```
python CRISPRtracrRNA.py --input_folder test_input_complete --anti_repeat_coverage_threshold 0.6
```

##### Ranking weights

Weights which can be used for the ranking of the results. The default values are:

| Description | Default value |
|-------------------------------------|-----|
| weight_crispr_array_score           | 0.5 |
| weight_anti_repeat_sim              | 0.5 |
| weight_anti_repeat_coverage         | 0.5 |
| weight_anti_sim_coverage            | 0.5 |
| weight_interaction_score            | 0.6 |
| weight_model_hit_score              | 0.9 |
| weight_terminator_hit_score         | 0.9 |
| weight_consistency_orientation      | 0.1 |
| weight_consistency_anti_repeat_tail | 0.1 |
| weight_consistency_tail_terminator  | 0.1 |
