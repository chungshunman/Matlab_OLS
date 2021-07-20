function output=my_summary(data)
    arguments 
        data table %A table of data set
    end
%%%%
% It generates a summary statistics table of the variables;
%
% INPUT data (required) (nobs x ncov) data set
%
% OUTPUT a summary table
% 
% USAGE  my_summary(data)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)
%%%%

colnames=data.Properties.VariableNames;
Mean=[]; Var=[]; Std=[]; Max=[]; Min=[]; %create empty variables
for j = 1:length(colnames)
   Mean(j)= mean(data{:, j}); %get the mean
   Var(j)=var(data{:, j}); %get the variance
   Std(j)=std(data{:, j}); %get the standard deviation
   Max(j)=max(data{:, j}); %get the maximum value
   Min(j)=min(data{:, j}); %get the minimum value
end

%put the results in a table
output=table(Mean', Var', Std', Max', Min');
output.Properties.VariableNames ={'Mean', 'Var', 'Std', 'Max', 'Min'};
output.Properties.RowNames=colnames;
disp('Summary Table')

