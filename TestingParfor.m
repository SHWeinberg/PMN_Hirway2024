Emat='test_100p_2501mcs_0tgfb_0.075spread1.25div_extraA0.96_extraBT0.6_ZEB1.125_SNAIL1_kDAE.scale1_base.mat'; % Epithelial Cells Confluence
gridmat='test_100p_1007mcs_0tgfb__ConfluenceK_kds_scale_0.5_BaseKupffer.mat'; % Using confluent K non transformed cells as base
loadE=load(Emat);
loadg=load(gridmat);
totalmat=cell(1,2);
totalmat{1}=gridmat;
totalmat{2}=Emat;

parval=cell(1,2);
rngcolumn=0;
tic
a=load(Emat);
for i=[1:5]
rngcolumn=rngcolumn+1;
parfor j=[1:5]

for k=1:3
    b=1+1;
end
end
end
toc