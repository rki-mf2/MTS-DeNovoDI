import os
import subprocess


directory = '/home/user/scratch/megaX/build/main/Meta_proteomics_Data/mgf/mgf_files/'
files = os.listdir(directory)

for file in files:
    mgf_file = directory + file
    command = "srun --gpus 1 --mem 50GB -c 32 python /home/user/scratch/SMSNet/run.py --model_dir /home/user/scratch/SMSNet/MassIVE_HCD --inference_input_file " + mgf_file + " --inference_output_file /home/user/scratch/SMSNet/output/"
    subprocess.run(command, shell=True)
