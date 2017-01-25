function S=com_dist(i,j)

persistent comp;

if isempty(comp)
    comp = textread('e_simmat_dc2.txt');
    comp = comp(:,2:(size(comp,2)));
end

S = 1.000000001-comp(i,j);