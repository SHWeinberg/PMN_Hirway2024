% 1-variable Hill curve
if 0
    E = 6;
    M = 1;
    phalf = 0.5; n = 6;
    
    pheno = linspace(0,1);
    A = E + (M-E)*pheno.^n./(pheno.^n + phalf^n);
    plot(pheno, A); grid on;
end

% 2-variable Hill curve
E = 1;
M = 5;
phalf = 0.5; n = 6;
pheno1 = linspace(0,1,50);
pheno2 = linspace(0,1,51);
[P1, P2] = meshgrid(pheno1, pheno2);

%% JcM
jcmE = 100; %OG was 2.5
jcmM = 10;%OG was 2.5
figure
n=6;
phalf=0.6;
J = jcmE + (jcmM-jcmE)*(P1).^n./((P1).^n + phalf^n);
plot(P1,J(1,:));

%% JcC
phenoparams.jcmE = 2.5; %OG was 2.5
phenoparams.jcmM = 2.5;%OG was 2.5
ratioE = 1; ratioM = 0.1;
pheno1 = linspace(0,1,50);
pheno2 = linspace(0,1,51);
[P1, P2] = meshgrid(pheno1, pheno2);
E = ratioE*phenoparams.jcmE;
M = ratioM*phenoparams.jcmM;
for i=3:3:9
    for phalf=0.4:.1:0.6
    n=i;
figure
% output governed by most mesenchymal-like
A_mostM = E + (M-E)*max(P1,P2).^n./(max(P1,P2).^n + phalf^n);
% subplot(1,2,1);
surf(P1, P2, A_mostM); title(['most mesenchymal, n:' num2str(n) ' phalf: ' num2str(phalf)])
xlabel('Cell P1'); ylabel('Cell P2'); zlabel('JCC');
    end
% union (similar to most mesenchymal-like)
% A_prod2 = E + (M-E)*(P1+P2-P1.*P2).^n./((P1+P2-P1.*P2).^n + phalf^n);
% subplot(1,2,2); surf(P1, P2, A_prod2); title('union');
% xlabel('Cell P1'); ylabel('Cell P2'); zlabel('JCC');
end

%%