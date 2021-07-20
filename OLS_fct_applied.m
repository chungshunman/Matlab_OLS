%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application of OLSmodel function
%%%%%%
% Authors: 
% Matthias Lukosch (20-601-050)
% Chung Shun Man (20-621-587)
%%%%%%
%% import data
[filename directory_name] = uigetfile('*.dat', 'Select a file');

fullname = fullfile(directory_name, filename);

data = readtable(fullname);

%% check whether there are missing values
sum(any(ismissing(data)),2)

%% Summary data

mysum=my_summary(data)


%% Data summary 1
%%%%%%%% a) regression
results = OLSmodel(data(:,1), data(:,{'MALE'}),...
    'h_robust', false);
results.fit_measures
%% Data summary 2
%%%%%%%% a) set up
FEMALE=1-table2array(data(:, 'MALE'));
data1=addvars(data,FEMALE);
%%%%%%%% b) regression
results = OLSmodel(data1(:, 'WAGE'),...
    data1(:, {'MALE'}));
%% Regression 1
%%%%%%%% a) set up
sqn=my_square(data, {'EXPER'}); data2=[data sqn];
%%%%%%%% b) regression
results = OLSmodel(data2(:, 'WAGE'), ...
                   data2(:, {'MALE', 'EDUC', 'EXPER', 'EXPER_2'})...
                   ,'resid_diagnostic', true);
%% Regression 2
%%%%%%%% a) set up
% transform variables
logn=my_log(data2,{'WAGE','EDUC', 'EXPER'}); data3=[data2 logn];
sqn=my_square(data3, {'EXPER_ln'}); data3=[data3 sqn];
data3.EXPER_ln=data3.EXPER_ln.*(data3.EXPER_ln~=-inf);
%%%%%%%% b) regression
results = OLSmodel(data3(:, 'WAGE_ln'), ...
          data3(:, {'MALE','EDUC_ln', 'EXPER_ln', 'EXPER_ln_2'}), ...
          'h_robust', true,'h_diagnostic', true);
%% Regression 3
%%%%%%%% a) set up
% create dummy variables
dummy=my_dummy(data, 'EDUC');
% create interaction terms
cross=my_cross(dummy, data(:,'EXPER'));
% create new data table
data4=[data(:,{'WAGE', 'MALE', 'EXPER'}),...
       dummy(:,2:end), cross(:,2:end)];
   
%%%%%%%%% b) regression
% apply OLSmodel
results=OLSmodel(data4(:, 'WAGE'), data4(:, {'MALE', 'EXPER', ...
    'EDUC2', 'EDUC3', 'EDUC4', 'EDUC5',...
    'EDUC2_EXPER', 'EDUC3_EXPER','EDUC4_EXPER', 'EDUC5_EXPER'}), ...
    'h_robust', true,'h_diagnostic', true, 'resid_diagnostic', true, 'spec_diagnostic', true);
