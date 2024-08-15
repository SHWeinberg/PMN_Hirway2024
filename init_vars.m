function def = init_vars(varargin)

if nargin>0
    for n=1:2:nargin
        switch varargin{n}
            case 'fmaflag'
                fmaflag=varargin{n+1};
            case 'amp'
                amp=varargin{n+1};
            case 'ratio'
                ratio=varargin{n+1};
            case 'divthresh'
                divthresh=varargin{n+1};
            case 'scurr'
                scurr=varargin{n+1};
            case 'ind'
                ind=varargin{n+1};
            case 'NVX'
                NVX = varargin{n+1};
            case 'NRINC'
                NRINC = varargin{n+1};
        end
    end
else
    fmaflag=0;
    ind=1;
    amp=1;
    ratio=2;
    scurr=rng('shuffle');
    divthresh=0.005;
    NVX = 100;
    NRINC = 1000; % 86400 for 3 days or 201,600 for 7 days
end

NULL  = 0;
FALSE  = 0;
TRUE = 1;

NVY = NVX; NV = NVX * NVY;
NNX = NVX + 1; NNY = NNX; NN  = NNX * NNY;
NDOF = 2 * NN;
VOXSIZE = 2.5E-6; % [m]
MAXNRITER = 1000;
ACCURACY = .00001;

% material properties
YOUNGS = 10E3; % [Pa]
POISSON = .45; %

% loading
LOAD = 0; % 3E6
FORCE = (LOAD*VOXSIZE);

% cells
IMMOTILITY = 1.0; %50
CELLDIAM = 2.0E-5; % cell diameter [m]
CELLRVOX = CELLDIAM/2/VOXSIZE; % cell radius [pixels]
TARGETVOLUME = 3.1415 * CELLRVOX * CELLRVOX; % targetvolume [pixels]
INELASTICITY = 500.0; % [-] 1.0E20 % [/m4]

NOSTICKJ = 500000.0; %10000 [/m] contact penalty for non-adhesive surface
JCM = amp * NOSTICKJ * VOXSIZE;  % cell-medium contact energy
JCC = ratio * JCM; % cell-cell contact energy

MAXDHSTR = 10.0; % unscaled at the moment,
THRESHOLDSTIFF = 15E3;    % threshold stiffness for durotaxis
STIFFSENSITIVITY = .0005; % steepness of durotaxis sigmoid
STIFFENINGSTIFF = .1; % steepness of strain-stiffening

THICKNESS = 10E-6; % effective thickness of substrate
CELLFORCE = 1.0E-5/THICKNESS; % [N]

SQ05 = .707107; %sqrt(.5), used often enough to make this convenient
RAND_MAX = 32767;

%binding kinetics
Ac = 2.5;
KOFF = 0.37;
Ac_Ka = 7.2E-7;
Ka = Ac_Ka / Ac;
KON = Ka * KOFF;
rhoR = 650; %molecules per um^2
ECAD_DENSITY = rhoR * Ac;
K2D = 3.1E-4; % 2D binding affinity um^2

dt = .01/max(KOFF, Ac*ECAD_DENSITY*ECAD_DENSITY*KON);

FIELDNAMES = {'fieldNames','fmaflag','ind','scurr','NULL','FALSE','TRUE','NVX','NVY','NV','NNX',...
    'NNY','NN','NDOF','VOXSIZE','NRINC','MAXNRITER','ACCURACY','YOUNGS','POISSON',...
    'LOAD','FORCE','IMMOTILITY','divthresh','CELLDIAM','CELLRVOX','TARGETVOLUME','INELASTICITY',...
    'NOSTICKJ','JCM','amp','JCC','ratio','MAXDHSTR','THRESHOLDSTIFF','STIFFSENSITIVITY',...
    'STIFFENINGSTIFF','THICKNESS','CELLFORCE','SQ05','RAND_MAX',...
    'ECAD_DENSITY','K2D','KON','KOFF','Ac','dt'};
def = v2struct(FIELDNAMES);

end
