function sqn=my_square(data, varnames)
    arguments 
        data table %A table of data set
        varnames string %The variable name
    end
%%%%
% It generates the second power of a variable;
%
% INPUT data (required) (nobs x ncov) data set
%
% OUTPUT sqn the square of a variable
% 
% USAGE  my_square(data, varnames)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)
%%%%

% Take the second power 
sqn= data{:, varnames}.^2;
sqn=array2table(sqn);
sqn.Properties.VariableNames= strcat(varnames,'_2');


end