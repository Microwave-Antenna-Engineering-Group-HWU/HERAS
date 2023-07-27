function property_array = extract_property(cell_array_of_objects,property_name)
% This function receives as input a cell array of objects belonging to
% different classes, which have a common property property_name. 
% The output property_array corresponds to an array containing the value
% associated to the property property_name for each element of the cell
% array.
% Author : Francesco Lisi
% Revisions : 0.1 - 03/03/2023

% Compute number of aelement in cell_array_of_objects
N_array=length(cell_array_of_objects);

% Compute the size of the value corresponding to the field property_name
eval(sprintf('property_size=size(cell_array_of_objects{1}.%s);',property_name));

% Extract dimension of the property
N_property_size=length(property_size);

% Extract all the cell_array_of_objects values corresponding to the field
% property_name in a matrix. The output matrix is obtained by stacking each
% cell array element field along the second dimension
eval(sprintf('property_array=cell2mat(cellfun(@(c) c.%s,cell_array_of_objects,''UniformOutput'',false));',property_name));

% shift dimension by 1 to the left
property_array=shiftdim(property_array,1);

% Extract property_size after the shiftdim operation
shifted_property_size=[property_size(2:end),property_size(1)];

% reshape and permute property_array so that the final size corresponds to
% [N_array,property_size]
property_array=permute(reshape(property_array,[shifted_property_size(1),N_array,shifted_property_size(2:end)]),[2,N_property_size+1,1,3:N_property_size]);

end