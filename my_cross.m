function cross=my_cross(dum, var)
    arguments 
        dum table %A table of dummy variables
        var table %A table of a covariate
    end
%%%%
% It creates cross terms with dummy variables and a covariate;
%
% INPUT dum (required) (nobs x ncov) dummy variables
%
%       var (required) (nobs x 1) a covariate
%
% OUTPUT the cross-terms
% 
% USAGE  cross=my_cross(dum, var)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)
%%%%

% Get the names of the dummy and the covariate
dum_name=dum.Properties.VariableNames;
var_name=var.Properties.VariableNames;

% Multiply both variables to create cross terms
cross=table2array(dum).*table2array(var);
cross=array2table(cross);
cross.Properties.VariableNames= strcat(string(dum_name), '_', var_name);

end