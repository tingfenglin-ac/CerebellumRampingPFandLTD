function cmap_new = interp_cmap(cmap_old, ncolor_new)
ncolor_old = size(cmap_old, 1);

if ncolor_old >= ncolor_new 
    cmap_new = cmap_old; 
    return;
end

pseu_vec_old = linspace(0, 1, ncolor_old); 

cmap_new = nan(ncolor_new, 3);
pseudo_vec_new = linspace(0, 1, ncolor_new); 

for i = 1:3
    cmap_new(:,i) = interp1(pseu_vec_old, cmap_old(:,i), pseudo_vec_new);
end

end