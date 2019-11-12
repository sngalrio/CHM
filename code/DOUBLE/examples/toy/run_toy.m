% run example of convex hull algorithm in toy model
% compute the projection in 2 dimensions - R1, R2

Aeq = load('toy_Aeq.txt');
beq = zeros(size(Aeq,1));

dims = [1 2]; % Indices of the reactions for the projection
dom = load('toy_domain.txt');
lbs = dom(:,1);
ubs = dom(:,2);

cd ../../chm
CH=computeCH(Aeq,lbs,ubs,dims);
