%% Write parameters to File
tic
excelfile='Test1.xlsx';
Collist=['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'];
Colnum=1;
Rownum=1;

%%
tic
SimID=cell(1,2);
SimID{1,1}='SimID';
SimID{1,2}=IDval;
othervalues=cell(10,2);
othervalues{1,1}='dim';
othervalues{1,2}=dimval;
othervalues{2,1}='beforeMCS';
othervalues{2,2}=beforemcsval;
othervalues{3,1}='afterMCS';
othervalues{3,2}=aftermcsval;
othervalues{4,1}='beforeTGFB';
othervalues{4,2}=beforetgfbval;
othervalues{5,1}='afterTGFB';
othervalues{5,2}=aftertgfbval;
othervalues{6,1}='beforedivideL';
othervalues{6,2}=pop_table.divL;
othervalues{7,1}='beforedivideH';
othervalues{7,2}=pop_table.divH;
othervalues{8,1}='afterdivideL';
othervalues{8,2}=pop_table.divL;
othervalues{9,1}='afterdivideH';
othervalues{9,2}=pop_table.divH;
othermodfn=fieldnames(othermod);
for i=1:numel(othermodfn)-4
    othervalues{i+9,1}=othermodfn{i};
    othervalues{i+9,2}=othermod.(othermodfn{i});
end
othervalues{i+10,1}=othermodfn{end};
othervalues{i+10,2}=othermod.(othermodfn{end});
writecell(SimID,excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum)]);
writecell(othervalues,excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+2)]);
%%
Colnum=1;
Rownum=24;
writematrix("TypeProb",excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum)]);
writematrix(prob_pop,excelfile,'Sheet',1,'Range',[Collist(Colnum+1) num2str(Rownum)]);
% writematrix("PopIndex",excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+1)]);
% writematrix(pop_index',excelfile,'Sheet',1,'Range',[Collist(Colnum+1) num2str(Rownum+1)]);
writematrix("CellParameters",excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+3)]);
writematrix(["Kupffer";"Stellate";"Tumor"],excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+4)]);
writematrix(celltypeparam,excelfile,'Sheet',1,'Range',[Collist(Colnum+1) num2str(Rownum+4)]);

%% 
Colnum=1;
Rownum=32;
Contactnames=cell(3,3);
Contactnames{1}='JcC';
Contactnames{2}=['K';'S';'T'];
Contactnames{3}=["K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T"];
Contactnames{4}='JcA';
Contactnames{5}=['K';'S';'T'];
Contactnames{6}=["K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T" "K" "S" "T"];
Contactnames{7}='JcM';
Contactnames{8}=['K';'S';'T'];
Contactnames{9}=["Less" "Mid" "More"];

writematrix(Contactnames{1},excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum)]);
writematrix(Contactnames{2},excelfile,'Sheet',1,'Range',[[Collist(Colnum) num2str(Rownum+1)] ':' [Collist(Colnum) num2str(Rownum+3)]]);
writematrix(Contactnames{3},excelfile,'Sheet',1,'Range',[[Collist(Colnum+1) num2str(Rownum)] ':' [Collist(Colnum) Collist(Colnum+1) num2str(Rownum+1)]]);
writecell(JccTab,excelfile,'Sheet',1,'Range',[Collist(Colnum+1) num2str(Rownum+1)]);
writematrix(Contactnames{4},excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+4)]);
writematrix(Contactnames{5},excelfile,'Sheet',1,'Range',[[Collist(Colnum) num2str(Rownum+5)] ':' [Collist(Colnum) num2str(Rownum+7)]]);
writematrix(Contactnames{6},excelfile,'Sheet',1,'Range',[[Collist(Colnum+1) num2str(Rownum+4)] ':' [Collist(Colnum) Collist(Colnum+1) num2str(Rownum+4)]]);
writecell(JcaTab,excelfile,'Sheet',1,'Range', [Collist(Colnum+1) num2str(Rownum+5)]);
writematrix(Contactnames{7},excelfile,'Sheet',1,'Range',[Collist(Colnum) num2str(Rownum+8)]);
writematrix(Contactnames{8},excelfile,'Sheet',1,'Range',[[Collist(Colnum) num2str(Rownum+9)] ':' [Collist(Colnum) num2str(Rownum+11)]]);
writematrix(Contactnames{9},excelfile,'Sheet',1,'Range',[[Collist(Colnum+1) num2str(Rownum+8)] ':' [Collist(Colnum+9) num2str(Rownum+8)]]);
writecell(JcmTab,excelfile,'Sheet',1,'Range',[Collist(Colnum+1) num2str(Rownum+9)]);

toc

%%
% getname(othermodfn{1})
% 
% toc
function s=getname(a)
s = inputname(1);
end


