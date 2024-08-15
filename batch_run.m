%LOAD PARAMETER FILE
load('sim_params3')
task = 1;
% task = str2num(getenv('SGE_TASK_ID'));

def=sim_params{task};

cpmfem('def',def)