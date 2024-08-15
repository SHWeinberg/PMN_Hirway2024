function D = material_matrix(D,def)
YOUNGS = def.YOUNGS;
POISSON = def.POISSON;

% D is pntr to stiffness matrix


planestress = true;

for m=1:3
    for n=1:3
        D(m,n) = 0;
    end
end

if (planestress)
    
    Es = YOUNGS/(1-POISSON*POISSON);
    % fill material matrix
    D(1+3*0) = Es * 1;
    D(2+3*1) = Es * 1;
    D(1+3*1) = Es * POISSON;
    D(2+3*0) = Es * POISSON;
    D(3+3*2) = Es * .5*(1-POISSON);
    
else % planestrain
    
    Es = YOUNGS/((1+POISSON)*(1-2*POISSON));
    % fill material matrix
    D(1+3*0) = Es * (1-POISSON);
    D(2+3*1) = Es * (1-POISSON);
    D(1+3*1) = Es * POISSON;
    D(2+3*0) = Es * POISSON;
    D(3+3*2) = Es * .5*(1-2*POISSON);
end
end
