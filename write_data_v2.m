function write_data_v2(filename, sim_data,def, cellmod, phenoparams, cell_data)
sim_data.def = def;
sim_data.cellmod = cellmod;
sim_data.phenoparams = phenoparams;
sim_data.cell_data = cell_data;
sim_data.name = filename;
save(filename,'-struct','sim_data');
gen_log(filename);

end
% function sim_file = store_sim(sim_data,def)
% jcm=def.amp*1E2;
% jcc=def.ratio*1E2;
% nrinc=def.NRINC;
% nvx=def.NVX;
% nvy=def.NVY;
% divthresh=def.divthresh*1E3;
% ind=def.ind;
% if def.fmaflag
%     fmaflag='atag';
% else
%     fmaflag='ctag';
% end
% 
% expression='sim_data_%s_jcm%d_jcc%d_mcs%d_%dx%d_div%d_rng%d';
% sim_file=sprintf(expression,...
%     fmaflag,jcm,jcc,nrinc,nvx,nvy,divthresh,ind);
% 
% sim_data.name=sim_file;
% sim_data.params=def;
% 
% end

function gen_log(sim_file)
longsec=toc;
simhour=floor(longsec/3600); simmin=floor(rem(longsec,3600)/60);
simsec=longsec-(3600*simhour+60*simmin);
simtime=sprintf('%.0fhr:%.fmin:%.f sec',simhour,floor(simmin),simsec);

meta_file=['metadata_', sim_file, '.txt'];

if ismac
    usrname=getenv('USER');
else
    usrname=getenv('username');
end

format_spec=[
    'user=%s\n'...
    'date=%s\n'...
    'MATLAB=%s %s\n'...
    'runtime=%s'...
    ];

fid=fopen(meta_file,'w');
fprintf(fid,format_spec,usrname,datetime,computer,version,simtime);
fclose(fid);
end

