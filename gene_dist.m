function S=gene_dist(i,j)

persistent target;

if isempty(target)
    target = textread('e_simmat_dg2.txt');
    target = target(:,2:(size(target,2)));
end

S = 1.000000001-target(i,j); 
