#!/bin/sh

### -- set the job Name --
#BSUB -J train
### -- ask for number of cores (default: 1) --
#BSUB -n 4
#BSUB -R "span[hosts=1]"
### -- specify queue -- voltash cabgpu gpuv100
#BSUB -q cabgpu
### -- set walltime limit: hh:mm --
#BSUB -W 120:00
### -- Select the resources: 1 gpu in exclusive process mode --:mode=exclusive_process
#BSUB -gpu "num=1:mode=exclusive_process"
## --- select a GPU with 32gb----
#BSUB -R "select[gpu40gb]"
### -- specify that we need 3GB of memory per core/slot --
#BSUB -R "rusage[mem=64GB]"
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o train.out
#BSUB -e train.err

# here follow the commands you want to execute

# submit with bsub < submit.sh
>train.out
>train.err

cd ~/projects/UVSNet/scripts
module load cuda/11.3
module load python3/3.7.11
module load pandas/1.3.1-python-3.7.11
module load numpy/1.21.1-python-3.7.11-openblas-0.3.17


pip install --upgrade pip
pip install  networkx
pip install --user matplotlib
pip install seaborn
pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113

#python3 /zhome/4c/8/164840/projects/UVSNet/scripts/01_simulate_compounds.py
#python3 /zhome/4c/8/164840/projects/UVSNet/scripts/02_simulate_strain_collection.py
#python3 /zhome/4c/8/164840/projects/UVSNet/scripts/04_train.py
python3 /zhome/4c/8/164840/projects/UVSNet/scripts/05_train_AL.py
python3 /zhome/4c/8/164840/projects/UVSNet/scripts/06_train_RL.py
python3 /zhome/4c/8/164840/projects/UVSNet/scripts/07_model_analysis.py