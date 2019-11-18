function cell_ind = cell_find(cell_contents, string)
% find the index location of string in cell_contents

cellfind = @(string)(@(cell_contents)(strcmp(string, cell_contents)));

cell_ind = find(cellfun(cellfind(string), cell_contents)==1);
