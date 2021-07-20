function dummies=my_dummy(data, var)
    arguments 
        data table %A table of the data set
        var string %The variable name
    end
%%%%
% It creates multiple dummy variables from an ordinal/discrete variable;
%
% INPUT data (required) (nobs x ncov) data set
%
% OUTPUT the dummy variables
% 
% USAGE  my_dummy(data, var)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)
%%%%

% get the unique values of the variable
uni=unique(data{:,var}); 
% create an empty variable to save the results
dummies=[]; 
[nobs ncov]=size(data);
% iterate through the observations to create dummy variables
for n = 1:nobs
    % iterate through the unique values of the variable
    for j = uni'
        if data{n, var}==j
            dummies(n,j)=1;
        else
            dummies(n,j)=0;
        end
    end
end

% turn the dummies into a table
dummies=array2table(dummies);
dummies.Properties.VariableNames=strcat(var, [string(uni)]');
end