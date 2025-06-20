srun --pty bash -i # goto interactive slurm session to avoid 10 minutes banhummer

eval __conda_setup="$('/usr/local/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"

conda activate depict

python2.7 /home/common/projects/depict_software/DEPICT_v194/src/python/depict.py for_depict_XXX.cfg
