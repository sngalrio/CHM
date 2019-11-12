% run example of convex hull algorithm in E.coli core model
% compute the projection in 2 dimensions - acetate and biomass

Aeq = load('ECC2_Aeq.txt');
beq = zeros(size(Aeq,1));

dims = [1 3]; % Biomass and acetate reaction indices
dom = load('ECC2_domain.txt');
lbs = dom(:,1);
ubs = dom(:,2);

cd ../../chm
CH=computeCH(Aeq,lbs,ubs,dims);
