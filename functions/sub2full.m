function fulldesign = sub2full(subdesign,dimfull)
% subdesign: subsystem design
% dimfull: dimensions of full system

% In the input data, say outx(:,:,1,1) 
% a column represents a y-direction material strip in FE model,
% a row represents a x-direction material strip in FE model.

sz = size(subdesign);
if length(sz)==2
    sz(3)=1;
end
dim_cell = dimfull./sz;

subdesign_cell = mat2cell(subdesign,ones(1,size(subdesign,1)),...
    ones(1,size(subdesign,2)),ones(1,size(subdesign,3)));

% do cell operation
fulldesign_cell = cellfun(@(x) x*ones(dim_cell),...
    subdesign_cell,'UniformOutput',false);

fulldesign = cell2mat(fulldesign_cell);

end

