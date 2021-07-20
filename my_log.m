function logn=my_log(data, varnames)
    arguments 
        data table %A table of a data set
        varnames string %The variable name
    end
%%%%
% It transforms an variable into its log form;
%
% INPUT data (required) (nobs x ncov) data set
%
% OUTPUT the variables in the log form
% 
% USAGE  my_log(data, varnames)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)
%%%%

% log the variable
logn=log(data{:, varnames}); 
logn=array2table(logn);
labels=strcat(varnames,'_ln');
logn.Properties.VariableNames=labels;

%%% Give a warning if there is infinity or minus infinity
[nobs ncov] = size(logn);

for v = 1:ncov
    ninf=sum(logn{:,v}==inf);
    nminf=sum(logn{:,v}==-inf);
    if nminf>0 | ninf>0
        code=sprintf('There is infinity or minus infinity in %s and is turned to zero',...
            string(labels(v)));
        disp(code);
    end
    logn{logn{:,v}==inf, v}=0;
    logn{logn{:,v}==-inf, v}=0;
end


end