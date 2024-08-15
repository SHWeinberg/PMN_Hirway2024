%LOAD PARAMETER FILE
load('sim_params4')
task = 251;
% task = str2num(getenv('SGE_TASK_ID'));

def=sim_params{task};

% change MCS
def.NRINC = 1000;

% change size
def.NVX = 200;
def.NVY = def.NVX; 
def.NV = def.NVX * def.NVY;
def.NNX = def.NVX + 1; 
def.NNY = def.NNX; 
def.NN  = def.NNX * def.NNY;
def.NDOF = 2 * def.NN;

cpmfem_v2(def)