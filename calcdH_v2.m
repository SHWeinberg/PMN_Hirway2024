function dH = calcdH_v2(ctag, ux, uy, csize, xt, xs, ...
    pick, ttag, stag, def, targetvols, jccs, jcms,phenotype,phenoparams,cellmod)

	dHcontact = 0;
	dHcontact = calcdHcontact_v3(ctag, xt, ttag, stag, def, jccs, jcms,phenotype,phenoparams,cellmod); %changed to v3 SUH 10/31

	dHvol = 0;
	dHvol = calcdHvol_v2(csize, ttag, stag, def, targetvols);

	dHstr = 0;
	dHstr = calcdHstrain(ux, uy, xt, xs, pick, ttag, stag, def);

	dH = dHcontact + dHvol + dHstr;

end
