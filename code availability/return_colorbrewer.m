function [cmap_mat, cmap_type] = return_colorbrewer(cmap_name, num_colors)
load('colorbrewer.mat', 'colorbrewer');
reverse_map = false;
reverse_suffix = '_r'; 

if strcmpi(cmap_name(end-length(reverse_suffix)+1:end), reverse_suffix)
    cmap_name(end-length(reverse_suffix)+1:end) = '';
    reverse_map = true;
end

available_cmaps = fieldnames(colorbrewer); 
cmap_ind = strcmpi(available_cmaps, cmap_name);

if ~any(cmap_ind)
    error('The map %s does not exist in the colorbrewer map, only %s', ...
        cmap_name, sprintf('"%s" ', available_cmaps{:}));
end

cmap_name = available_cmaps{cmap_ind};
cmap_struct = colorbrewer.(cmap_name);

cmap_type = cmap_struct.type; 
cmap_numvec = cmap_struct.num; 
max_cnum = max(cmap_numvec);
min_cnum = min(cmap_numvec); 

if ~exist('num_colors', 'var') 
    num_colors = max_cnum; 
end

if num_colors < min_cnum
    error('The colorbrewer map needs to have at least %d colors', min_cnum);
end


if num_colors <= max_cnum
    cmap_mat = cmap_struct.(sprintf('x%d',num_colors));
else
    cmap_mat = cmap_struct.(sprintf('x%d',max_cnum));
    if strcmpi(cmap_type, 'qual') 
        num_repeat = ceil(num_colors/max_cnum); 
        cmap_mat = repmat(cmap_mat, [num_repeat, 1]); 
        cmap_mat = cmap_mat(1:num_colors,:); 
    else
        cmap_mat = interp_cmap(cmap_mat, num_colors);
    end 
end

if reverse_map
    cmap_mat = flipud(cmap_mat); 
end