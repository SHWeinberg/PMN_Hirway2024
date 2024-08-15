load('PMNrandom.mat')
xrange=[1:3000]*4.8/60;
figure
subplot(1,3,1)
plot(xrange,PMNrandom{1,2}, 'LineWidth', 2);
hold on
plot(xrange,PMNrandom{2,2},'LineWidth',2);
hold on
plot(xrange,PMNrandom{3,2}, 'LineWidth', 2);
hold on
plot(xrange,PMNrandom{4,2}, 'LineWidth', 2);

xlabel('Time (hours)','FontSize',24);
ylabel('ECM-Bound TGF-?  (?M)','FontSize',24);
title('Average ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;

subplot(1,3,2)
plot(xrange,(PMNrandom{1,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot(xrange,(PMNrandom{2,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot(xrange,(PMNrandom{3,3}/(100*100))*100, 'LineWidth', 2);
hold on
plot(xrange,(PMNrandom{4,3}/(100*100))*100, 'LineWidth', 2);

legend(PMNrandom{1,1},PMNrandom{2,1},PMNrandom{3,1},PMNrandom{4,1},'FontSize',18);
xlabel('Time (hours)','FontSize',24);
ylabel('Grid Percentage (%)','FontSize',24);
title('90% Max ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;

subplot(1,3,3)
plot([0.1, 1, 5, 10],[(PMNrandom{1,3}(end)/(100*100))*100,(PMNrandom{2,3}(end)/(100*100))*100,(PMNrandom{3,3}(end)/(100*100))*100,(PMNrandom{4,3}(end)/(100*100))*100],'k-o', 'LineWidth', 2);
xlabel('TGF-? Dose (?M)','FontSize',24);
%ylabel('ECM-Bound TGF-?  (?M)','FontSize',24);
ylabel('Grid Percentage (%)','FontSize',24);
title('Final 90% Max ECM-Bound TGF-?','FontSize',26);
ax = gca;
ax.FontSize = 18;
set(gca, 'XScale', 'log')
