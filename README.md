# PHH-metaPRS
PHH-metaPRS (Predicting Human Height with meta-PRS) is  a meta-PRS height prediction model for men and women, which included nine sets of PRSs, namely three sets of PRSs directly calculated by three GIANT Height PRS weight (‘ALL’, ‘EAS’, ‘EUR’), three sets of PRSs calculated by PRScs based on three GIANT Height GWAS (‘ALL’, ‘EAS’ for men, ‘EUR’, ‘HIS’ for women), two sets of PRSs calculated by PCNN based on a GIANT Height GWAS (‘ALL’ for men, ‘EAS’ for women) and a Biobank Japan Height GWAS (‘BBJ’), and one set of PRS calculated by PLR based on a Wegene height GWAS. 

Note: The model only works when the coverage between variants and summary statistics is larger than 0.4. If you want the prediction value even if the coverage rate is low, you can modify the default threshold value of ```coverage_threshold``` in 'height_pred.r'.



## Installation:
### From source
Download a local copy of :
```
git clone https://github.com/zzfcharlie/PHH-metaPRS.git
```
### Dependencies:
Python: Pytorch-cpu, Scikit-learn, Numpy, Pandas. 

R: bigsnpr, dplyr, data.table, Matrix, doParallel, recticulate, and all of their dependencies.

## An example to make prediction on testing data:

### STEP 1: Set up Python environment and install R packages.
We recommend you use conda to manage dependencies in different environments. If Conda hasn't been installed on your system yet, please visit https://www.anaconda.com/download for detailed installation information. Our analyses are currently conducted on CPU, and we are considering using GPU devices to train PyTorch model in the future. Refer to https://pytorch.org/get-started/locally/ for instructions on installing PyTorch, with steps varying based on your operating system.

First, create and activate a conda environment and install the Python packages. (Please skip this step if you've already done this before.)
```
# The example is performed under Windows 10.
conda create -n yourenv python=3.9
conda activate yourenv
conda install pytorch torchvision torchaudio cpuonly -c pytorch
conda install scikit-learn pandas numpy
```
Please record your conda environment name. We use ```yourenv``` in our example.

Then, install R packages as follows:
```R
install.packages('bigsnpr')
install.packages('dplyr')
install.packages('doParallel')
install.packages('recticulate')
install.packages('data.table')
```


### STEP 2: Load dataset (PLINK files or vcf file).
```R
setwd('PHH-metaPRS')
source('read_file.r')
source('height_pred.r')
#Our example data is in PLINK format.
obj.bigSNP <- read_file(
  input_file_dir = 'data/testing',
  #plink_dir = 'plink',
  #bfile_out_dir = 'path/to/store_bfile'
)
G <- obj.bigSNP$genotypes
map <- obj.bigSNP$map[-3]
names(map) <- c('chr', 'rsid', 'pos', 'a1','a0')

```
* ```read_file()``` accepts two types of input: VCF file(.vcf) and PLINK files (.bed, .bim, .fam) without extension.
* ```bfile_out_dir``` refers to the path where to store PLINK files and is only enabled when the input is in VCF format. Moreover, if you don't specify ```bfile_out_dir```, PLINK files will be generated under the same directory as your input file. 
* ```G``` refers to genotyped variants coded in '0/1/2'.
* ```map``` refers to variants information.
  
| chr | rsid        | pos   | a1 | a0 |
| --- | ----------- | ----- | -- | -- |
| 1   | rs13303291  | 862093| T  | C  |
| 1   | rs4040604   | 863124| G  | T  |
| ... | ...         | ...   | ...| ...|
| 22  | rs5771014   | 51216731 | C  | T  |
| 22  | rs28729663  | 51219006 | A  | G  |
| 22  | rs9616978   | 51220319 | G  | C  |

where chr, pos, a0, and a1 represent chromosome, position, reference allele, and alternative allele respectively.
* Our model doesn't handle missing values. You can use ```snp_fastImpute()``` or ```snp_fastImputeSimple()``` in ```bigsnpr``` to impute missing values of genotyped variants before the prediction.

### STEP 3: Predicting Height using Meta-PRS model.
```R
result <- height_predict(
  G,
  map,
  gender_ID_list_dir = 'data/gender_ID_list.txt',
  match_by_pos = TRUE,
  height_material_dir = 'path/to/height_material',
  env_name = 'yourenv',
  result_dir = NULL,
  Ncores = 2
)
```
* ```G``` refers to genotyped variants coded in '0/1/2'.
* ```map``` refers to variants information.
* ```gender_ID_list_dir``` : The path of a text file recording the individual's ID and gender(1 represents men and 2 represents women), where the file looks like:

| ID                                | gender |
|-----------------------------------|--------|
| 000861b9d9b502b5939136c8a8e66d0b | 2      |
| 001ca04596d464db83a8d3747a3167be | 2      |
| 002d43d7e8a42d32b8bd668e3dbc40b2 | 2      |
| 003e3c8006e11f9bb2faf166a336f898 | 1      |
| 004e0559b7d68acf4fdd7fa44aa31778 | 1      |
| 007da44b9de4ad3931efdf647b938e02 | 2      |
| 009fdd5194b699da3be445d4714f5860 | 2      |
| 009fe3c2d4b3ca29698ed595690c6ccc | 1      |
| 0012c8fab85c015c975bed170505eb3a | 2      |
| 0014c58204d4c306a8d3019eafdf4a57 | 1      |


* ```match_by_pos``` : If TRUE, variants were matched by position; otherwise, variants were matched by rsid. Default is TRUE.
* ```height_material_dir``` : 'height_material' contains the models required for the prediction process. Please download the 'height_material' from 'www.abc.com' to your local device and record the directory path.
* ```env_name``` : The name of the conda environment we created earlier, such as 'yourenv' in our example.
* ```result_dir``` : The path to save the result. If you don't specify, the result will be generated under your current R workspace.
* ```Ncores``` represents the number of CPU cores used for parallel processing.

## Final result:
The prediction result would be stored in 'result_dir/height_result.txt' with 3 columns (ID, gender, and the height prediction value), which looks like:
| ID                                | gender | height |
|-----------------------------------|--------|--------|
| 000861b9d9b502b5939136c8a8e66d0b | 2      | 159.1  |
| 001ca04596d464db83a8d3747a3167be | 2      | 167.8  |
| 002d43d7e8a42d32b8bd668e3dbc40b2 | 2      | 163.5  |
| 003e3c8006e11f9bb2faf166a336f898 | 1      | 180.2  |
| 004e0559b7d68acf4fdd7fa44aa31778 | 1      | 177.5  |
| 007da44b9de4ad3931efdf647b938e02 | 2      | 160.3  |
| 009fdd5194b699da3be445d4714f5860 | 2      | 161.2  |
| 009fe3c2d4b3ca29698ed595690c6ccc | 1      | 175.8  |
| 0012c8fab85c015c975bed170505eb3a | 2      | 163.6  |
| 0014c58204d4c306a8d3019eafdf4a57 | 1      | 177.5  |

