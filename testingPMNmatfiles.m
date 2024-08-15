load('test_100p_800mcs_10tgfb__BeforeTumor_divf_0.5_param_14_scale_0.5ECMscaled0_1000937.mat');
%%
figure
avgp=zeros(800,1);
for i=1:800
    avgp(i)=mean(phenotype{i});
end
plot(avgp)