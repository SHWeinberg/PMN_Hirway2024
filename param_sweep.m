%parameter vectors
% first analysis 
fmaflag=[0 1];                      %FMA calculation single (0) multicell (1)
% amp=[1 2 10];                       %Amplify Jcm and Jcc
% ratio=[0.5 1 2 10];                 %Jcc:Jcm
divthresh=[0.005 0.002];   %minimum division probability
% Nsamples=5;  %sample size

amp=[1 2 3 ];                       %Amplify Jcm and Jcc
ratio=[0.5 1 1.5 2 2.5 3 ];                 %Jcc:Jcm
Nsamples=5;      

%default simulation parameters
NVX = 100;
NRINC = 1000;

La=length(amp);
Lr=length(ratio);
Ld=length(divthresh);

%rng set
scurr=cell(Nsamples,1);
for xi=1:Nsamples
    scurr{xi}=rng('shuffle');       %Random seed for each sample
    pause(0.1)
end

%parameter set
sim_params=cell(2*Nsamples*La*Lr*Ld,1);
n=0;
for xf=1:2
    for xi=1:La
        for xj=1:Lr
            for xd=1:Ld
                for xn=1:Nsamples
                    n=n+1;
                    sim_params{n}=init_vars(...
                        'fmaflag',fmaflag(xf),...
                        'amp',amp(xi),...
                        'ratio',ratio(xj),...
                        'divthresh',divthresh(xd),...
                        'scurr',scurr{xn},...
                        'ind',xn,...
                        'NVX',NVX,...
                        'NRINC',NRINC...
                        );
                end
            end
        end
    end
end

save('sim_params4','sim_params')
