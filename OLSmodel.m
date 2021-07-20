function ols = OLSmodel(y_table,cov_table,options)
    arguments 
        y_table table 
        cov_table table     
        options.h_robust logical {mustBeNumericOrLogical} = true
        options.h_diagnostic logical {mustBeNumericOrLogical} = false
        options.resid_diagnostic logical {mustBeNumericOrLogical} = false
        options.spec_diagnostic logical {mustBeNumericOrLogical} = false
    end
% Description OLSmodel
% OLS regression and diagnostics of OLS assumptions.
% 
% INPUT y_table (required)  (nobs x 1)  dependent variable
%
%       cov_table (required)  (nobs x ncov) independent variables
%
%       h_robust  (optional) logical  Specify whether to use
%                                     heteroskedasticity robust standard errors
%                                     (default = true)
%
%       h_diagnostic (optional) logical Specify whether to test for
%                                       heteroskedasticity (default = false)
%       
%       resid_diagnoistic (optional) logical Specify whether to perform
%                                            further OLS residual diagnostics.
%                                            (default = false)
%       spec_diagnostic (optional) logical Specify whether to perform a
%                                          Regression Specification Error Test (RESET).
%                                          (default = false)
% OUTPUT structure array (1x1)
% 
% USAGE  OLSmodel(y,cov, 'h_robust', false)
%
% Authors: Matthias Lukosch (20-601-050)
%          Chung Shun Man (20-621-587)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1. Nested necessary functions (could also be stored externally)
% cumulative distribution function student's t distribution
function cdfs = t_cdf(t,df)

% define pdf
    function y = t_pdf2(x)
        y = (exp(gammaln((df+1)./2) - gammaln(df./2))) ./ ...
            (sqrt(df .* pi)*(1+(x.^2)./df).^((df+1)./2));
    end

% determine count of t statistics
ct = size(t,1);

% preallocate storage
cdfs = [zeros(ct,1)];

% loop over the values of the t statistic
for k = 1:numel(t)
    cdfs(k) = integral(@t_pdf2, -Inf, t(k), 'RelTol', 1e-40);
end

% correct precision
cdfs(cdfs > 1) = 1;

end

%%%% cumulative distribution function F distribution
function cdf = f_cdf(t, df1, df2)

% define pdf
    function y = f_pdf2(x)
        z = gammaln((df1 + df2)./2) - (gammaln(df1./2)+gammaln(df2./2))...
            + df1./2 .* (log(df1)-log(df2));
        k = ((df1-2)./2) .* log(x) - ((df1 + df2)./2) .*...
            log(1+(df1./df2).*x);
        y = exp(z + k);
    end

% compute integral
cdf = integral(@f_pdf2, 0, t, 'RelTol', 1e-40);

% correct precision
cdf(cdf > 1) = 1;

end

%%%% cumulative distribution function chisquare distribution
function cdf = chi2_cdf(t,df)
    
    % define pdf
    function y = chi2_pdf(x)
        a = ((df-2)./2) .* log(x) - x./2;
        b = (df./2) .* log(2) + gammaln(df./2);
        y = exp(a-b);
    end
    
    % compute integral
    cdf = integral(@chi2_pdf, 0, t);
    
    % correct precision
    cdf(cdf > 1) = 1;    
end

%%%% basic input check
function check_input_basic(nobs, nobs2, ncov, cov_table, y_table)
    % # observations of both inputs        
    if ~isequal(nobs, nobs2)
            ME_1 = MException('OLSmodel:inputError', 'The dependent variable must have as many observations as the covariate(s).');
            throwAsCaller(ME_1)
    end
    % # observations vs. # covariates
    if nobs < ncov
            ME_2 = MException('OLSmodel:inputError', 'There are less observations than covariates.');
            throwAsCaller(ME_2)
    end
    % check whether there are missing values in the table of covariates
    if sum(any(ismissing(cov_table)),2) ~= 0
            ME_3 = MException('OLSmodel:inputError', 'The covariate table contains missing values.');
            throwAsCaller(ME_3)
    end
    if any(ismissing(y_table)) ~= 0
            ME_4 = MException('OLSmodel:inputError', 'The dependent variable table contains missing values.');
            throwAsCaller(ME_4) 
    end
end

%%%% detect perfect multicollinearity
function check_perfect_multicollinearity(x)
    % x'x is not invertible if x'x is singular
    if rank(x) < ncov + 1
            ME = MException('OLSmodel:inputError', 'Perfect multicollinearity detected');
            throwAsCaller(ME)
    end
end

%%%% evaluate imperfect multicollinearity
function corr = corr_matrix(x,ncov)
    corr=[zeros(ncov)];
    for n=2:(ncov+1) %row
        for j=2:(ncov+1) %column
            num= sum((x(:,n)-mean(x(:,n))).*(x(:,j)-mean(x(:,j)))); %numerator
            den1=sqrt(sum((x(:,n)-mean(x(:,n))).^2)); %first denominator
            den2=sqrt(sum((x(:,j)-mean(x(:,j))).^2)); %second denominator
            corr((n-1),(j-1))=num/(den1*den2); %get the Pearson correlation
            if n>j % used to avoid double warnings
                if corr((n-1),(j-1))==-0.7 | corr((n-1),(j-1))==0.7 && n~=j
                    warn='There is perfect multicollinearity between %s and %s';
                    E=sprintf(warn, string(labels(n)), string(labels(j)));
                    disp(E) %display the warning
                end
            end
        end
    end
end

%%%% visual inspection heteroskedasticity
function h_visual_inspection(ncov, x, label_var, resid, y_head)
    if ncov == 1
        figure;
        scatter(x(:,2), resid, 'filled');
        xlabel(label_var(1));
        ylabel('Residuals');
        title('Heteroskedasticity Inspection');
    elseif ncov > 1
        figure;
        scatter(y_head, resid, 'filled');
        xlabel('Fitted Values');
        ylabel('Residuals');
        title('Heteroskedasticity Inspection');
    end
end

%%%% histogram of OLS residuals
function my_hist(resid)
        figure;
        % create histogram of residuals
        histogram(resid);
        % add x axis label
        xlabel('OLS residuals');
        % add y axis label
        ylabel('Frequency');
        % add title
        title('Histogram of OLS residuals');      
end

%%%% scatterplot of OLS residuals vs. lagged OLS residuals
function my_autocorr(resid)
     % create vector with lag 1 observations
     lagged_resid = [NaN; ols.resid];
     lagged_resid(end) = [];
     
     % plot OLS residuals vs. OLS lagged residuals
     figure;
     scatter(lagged_resid, resid, 'filled');
     xlabel('Lagged OLS residuals (l = 1)');
     ylabel('OLS residuals');
     title('Autocorrelation Inspection');    
end

%%%% Breusch-Pagan-Godfrey Test (Heteroskedasticity)
function BPG = breusch_pagan_godfrey_test(resid, x, nobs, ncov)
    % dependent variable: squared OLS residuals
    BPG.resid_sq = resid.^2;
    
    % auxiliary regression coefficients
    BPG.BPG_beta = (x' * x)\(x' * BPG.resid_sq);
    
    % predicted squared residuals
    BPG.resid_sq_head = x * BPG.BPG_beta;
    
    % coefficient of determination R^2
    BPG.BPG_r_2 = (BPG.resid_sq'*x*BPG.BPG_beta - nobs * mean(BPG.resid_sq)^2)/ (BPG.resid_sq'*BPG.resid_sq - nobs * mean(BPG.resid_sq)^2);
    
    % calculate test statistic
    BPG.BPG_statistic = nobs * BPG.BPG_r_2;
    
    % p-value
    BPG.p_value_BPG = 1 - chi2_cdf(BPG.BPG_statistic, ncov);
    
    % print result at 5 % significance level
    if BPG.p_value_BPG <= 0.05
        disp(' ');
        disp('Breusch–Pagan–Godfrey Test:');
        test=sprintf('Test statsitics: %d', BPG.BPG_statistic);
        disp(test);
        pv=sprintf('P-value: %d', BPG.p_value_BPG);
        disp(pv);
        disp('The null hypothesis of homoskedasticity can be rejected at the 5 % significance level.');
        disp(' ');
    else 
        disp(' ');
        disp('Breusch–Pagan–Godfrey Test:');
        test=sprintf('Test statsitics: %d', BPG.BPG_statistic);
        disp(test);
        pv=sprintf('P-value: %d', BPG.p_value_BPG);
        disp(pv);
        disp('The null hypothesis of homoskedasticity cannot be rejected at the 5 % significance level.');
        disp(' ');
    end
end

%%%% Jarque-Bera Test (Distribution OLS residuals)
function JB = jarque_bera_test(resid, nobs)
     % estimator of symmetry
     JB.sym_head = (mean(resid.^3)) / ((mean(resid.^2))^(3/2));
     % estimator of kurtosis
     JB.kur_head = (mean(resid.^4)) / ((mean(resid.^2))^2);
     
     % JB test statisic
     JB.JB_statistic = nobs * (((JB.sym_head)^2) / 6 + ((JB.kur_head - 3)^2)/ 24);
     
     % p- value
     JB.p_value_JB = 1 - chi2_cdf(JB.JB_statistic,2);
     
     % print results
     if JB.p_value_JB <= 0.05
         disp('Jarque-Bera Test:');
         disp('The null hypothesis of normally distributed residuals can be rejected at the 5 % significance level.');
     else
         disp('Jarque-Bera Test:');
         disp('The null hypothesis of normally distributed residuals cannot be rejected at the 5 % significance level.');
         disp(' ');
     end
end

%%%% Durbin-Watson Autocorrelation Test Statistic
function DW = durbin_watson(resid)
     % create vector with lag 1 observations
     lagged_resid = [NaN; resid];
     lagged_resid(end) = [];
     lagged_resid(1) = [];
     
     % create residuals vector starting at N = 2
     resid_N2 = resid(2:end);
     
     % estimator of autocorrelation in AR(1) process
     DW.cor_head = (lagged_resid' * resid_N2) / (sum(lagged_resid.^2));
     
     % test statistic
     DW.DW_statistic = (sum((resid_N2 - lagged_resid).^2)) / (sum(resid.^2));
     
     % print results
     disp('Please use the DW statistic that is given below to perform DW tests.')
     disp(['Durbin-Watson Test statistic: ', num2str(DW.DW_statistic)]);
     disp(['Estimate autocorrelation (l=1): ', num2str(DW.cor_head)]);
     disp(' ');     
end

%%%% RESET (Regression Specification)
function reset = RESET(x, y, y_head, ncov, nobs)
    % set covariates
    reset.x = [x, y_head.^2, y_head.^3, y_head.^4];
    % set # covariates
    reset.ncov = ncov + 3;
    % reset dgf
    reset.dgf = nobs - reset.ncov - 1;
    % calculate OLS estimates
    reset.beta = (reset.x' * reset.x)\(reset.x' * y);
    % calculate predicted values
    reset.y_head = reset.x * reset.beta;
    % calculate residuals
    reset.resid = y - reset.y_head;
   
    % calculate residual variance
    reset.resid_variance = (reset.resid' * reset.resid) / reset.dgf;
   
    
    % F test that all coefficients on higher orders are statistically not
    % different from zero
    % F-test statistic
    R = [zeros(3, ncov + 1), eye(3,3)];
    q = zeros(3,1);
    
    reset.F_statistic = (((R*reset.beta - q)' * inv(R*inv(reset.x'*reset.x)*R')*...
    (R*reset.beta - q)) / reset.ncov) / (reset.resid_variance);

    % p - value of F statistic
    reset.p_valueF = 1 - f_cdf(reset.F_statistic, reset.ncov, reset.dgf);
    
    % print results
    if reset.p_valueF <= 0.05
        disp('The null hypothesis of a linear specification can be rejected at the 5 % significance level.');
        disp(' ');
    else
        disp('The null hypothesis of a linear specification cannot be rejected at the 5 % significance level.');
        disp(' ');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  2. Read input and first checks

% get the outcome variable label
label_y = y_table.Properties.VariableNames;
% get the covariate values in matrix form
y = y_table{:,label_y};
% get the covariate labels
label_var =cov_table.Properties.VariableNames;
% get the covariate values in matrix form
cov = cov_table{:,label_var};
% add the label for the intercept
labels =['Intercept', label_var];

% store # observations and # covariates
[nobs ncov] = size(cov_table); nobs2 = size(y,1);
ols.nobs = nobs;
ols.ncov = ncov;

% get number of arguments passed
nargin;

% check inputs
check_input_basic(nobs, nobs2, ncov, cov_table, y_table);

% degress of freedom
dgf = nobs - ncov - 1;

% add constant
x = [ones(nobs,1), cov];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3. Multicollinearity Check

% Multicollinearity
% perfect
% check whether covariates are linearly dependent (perfect multicollinearity)
check_perfect_multicollinearity(x);

% imperfect
% check pearson correlation among regressors (only a sufficient condition for multicollinearity)
ols.cov_corr = array2table(corrcoef(cov));
ols.cov_corr.Properties.VariableNames = label_var;
ols.cov_corr.Properties.RowNames = label_var;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 4. OLS implementation
% OLS estimator
ols.beta = (x' * x)\(x' * y);

% OLS prediction
ols.y_head = x * ols.beta;

% OLS residuals
ols.resid = y - ols.y_head;

% estimator of the variance of the error term
ols.resid_variance = (ols.resid' * ols.resid)/dgf;

% homoskedasticity vs. heteroskedasticity
% homoskedasticity
if options.h_robust == false
    % variance-covariance matrix of the OLS estimator
    ols.homo.covariance_beta = ols.resid_variance * inv(x'*x);
    % variance OLS estimator
    ols.homo.variances_beta = diag(ols.homo.covariance_beta);
    % standard error OLS estimator
    ols.homo.se_beta= sqrt(ols.homo.variances_beta);
    
    % hypothesis testing
    % t-test statistic
    ols.homo.t_statistic = ols.beta ./ ols.homo.se_beta;

    % p - values
    % pdf of the Student's t distribution
    ols.homo.p_valueT = (1 - t_cdf(abs(ols.homo.t_statistic),dgf)) .* 2;

% heteroskedasticity
elseif options.h_robust == true
    
    % Huber-White sandwich estimator
    % non-constant variance
    ols.hetero.resid_omega=ols.resid.^2/ols.resid_variance;
    % heteroskedasticity robust variances
    ols.hetero.variance_beta = diag(ols.resid_variance *...
        inv(x'*x)*x'*diag(ols.hetero.resid_omega)*x*inv(x'*x));
    % heteroskedasticity robust standard errors
    ols.hetero.se_beta=sqrt(ols.hetero.variance_beta);    

    % hypothesis testing
    % t-test statistic
    ols.hetero.t_statistic = ols.beta ./ ols.hetero.se_beta;
    
    % p-values
    ols.hetero.p_valueT = (1 - t_cdf(abs(ols.hetero.t_statistic),dgf)).* 2;
end

% F-test statistic
R = [zeros(ncov,1), eye(ncov)];
q = zeros(ncov,1);
ols.F_statistic = (((R*ols.beta - q)' * inv(R*inv(x'*x)*R')*...
    (R*ols.beta - q)) / ncov) / (ols.resid_variance);

% p - value of F statistic
ols.p_valueF = 1 - f_cdf(ols.F_statistic, ncov, dgf);

% measures of fit
% sum of squared residuals
ols.fit_measures.sum_squared_resid = ols.resid' * ols.resid;
% standard error of regression
ols.fit_measures.se_reg = sqrt(ols.resid_variance);

% coefficients of determination
% multiple R squared
ols.fit_measures.r_2 = (y'*x*ols.beta - nobs * mean(y)^2)/ (y'*y - nobs * mean(y)^2);
% adjusted R squared
ols.fit_measures.r_2_adjusted = 1-(1-ols.fit_measures.r_2)* nobs/dgf;
% mean squared error
ols.fit_measures.mse = mean((y - ols.y_head).^2);
% mean absolute error
ols.fit_measures.mae = mean(abs(y - ols.y_head));

% information criteria
% Akaike's information criterion
ols.info_criteria.aic = log(ols.fit_measures.sum_squared_resid/nobs) + (2*(ncov+1))/nobs;
% Schwarz criterion
ols.info_criteria.schwarz = log(ols.fit_measures.sum_squared_resid/nobs) + (ncov * log(nobs))/nobs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 5. Diagnostics

%%% Heteroskedasticity Diagnostic
if options.h_diagnostic == true

    % describe output
    disp('################################################################################');
    disp('################################################################################');
    disp('Heteroskedasticity Diagnostic:');
    disp('Visual inspection: Please see the created figure.');
    disp('Formal Test: Breusch-Pagan-Godfrey Test.');
    
    % visual inspection
    h_visual_inspection(ncov, x, label_var, ols.resid, ols.y_head);
    
    % formal test
    % Breusch–Pagan–Godfrey Test
    ols.diagnostics.BPG = breusch_pagan_godfrey_test(ols.resid, x, nobs, ncov);
    
end    
    
%%%% further OLS residual diagnostics
if options.resid_diagnostic == true
     
     % describe output
     disp('################################################################################');
     disp(' ');
     disp('Further OLS Residuals Diagnostics:');
     disp('Visual inspection: Please see the created figure.');
     disp('Formal Tests:');
     disp(' ');
     
     %%%% Distribution of Residuals
     % visual inspection
     % create histogram of OLS residuals
     my_hist(ols.resid);
     
     % print info
     disp('Regarding distribution of OLS residuals:');
     % formal Jarque-Bera Test
     ols.diagnostics.JB = jarque_bera_test(ols.resid, nobs);
     
     %%%% Autocorrelation
     % visual inspection
     my_autocorr(ols.resid);
     
     % print info
     disp(' ');
     disp('Regarding autocorrelation (lag = 1):');
     % formal Durbin-Watson Test statistic
     ols.diagnostics.DW = durbin_watson(ols.resid);     
end

%%%% Functional specification diagnostics

if options.spec_diagnostic == true
    
    disp('################################################################################');
    disp(' ');
    disp('Regression Specification Error Test (RESET): ');
    % Regression specification error test (RESET)
    ols.diagnostics.reset = RESET(x, y, ols.y_head, ncov, nobs);
       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 6. Output
disp('################################################################################');
disp(' ');
disp(['Execution Date:', string(datetime('now'))]);
disp('################################################################################');
disp(' ');
disp('Regression Summary Table:');
disp(' ');
Dep=sprintf('Dependent variable:  %s', string(label_y));
disp(Dep);
disp(' ');

% heteroskedasticity output table
if options.h_robust == true
    % combine column vectors with results
    ols.hetero.table = round([ols.beta, ols.hetero.se_beta,...
        ols.hetero.t_statistic, ols.hetero.p_valueT],6);
    % convert array to table
    ols.hetero.table = array2table(ols.hetero.table);
    % set column names
    ols.hetero.table.Properties.VariableNames = ["Coefficient",...
        "SE", 't-Statistic', 'p-Value'];
    % set row names
    ols.hetero.table.Properties.RowNames = labels;
    % display output table
    disp(ols.hetero.table);

% homoskedasticity output table
elseif options.h_robust == false
    % combine column vectors with results
    ols.homo.table = round([ols.beta, ols.homo.se_beta, ...
        ols.homo.t_statistic, ols.homo.p_valueT],6);
    % convert array to table
    ols.homo.table = array2table(ols.homo.table);
    % set column names
    ols.homo.table.Properties.VariableNames = ["Coefficient",...
        "SE", 't-Statistic', 'p-Value'];
    % set row names
    ols.homo.table.Properties.RowNames = labels;
    % display output table
    disp(ols.homo.table);
end 


% output table measures of fit
disp('################################################################################');
disp(' ');
disp(' ');
disp('Measures of Fit:');

% array
measures_of_fit = round([ols.fit_measures.sum_squared_resid; ols.fit_measures.se_reg;
                   ols.fit_measures.r_2; ols.fit_measures.r_2_adjusted;
                   ols.fit_measures.mse; ols.fit_measures.mae],6);
% convert array to table             
measures_of_fit = array2table(measures_of_fit);

% set row names
measures_of_fit.Properties.RowNames = {'Sum of squared residuals', ...
    'SE of regression','R-squared', 'Adjusted R-squared', 'MSE', 'MAE'};
                                    
% set column name
measures_of_fit.Properties.VariableNames = {' '};
                             
% print results                                   
disp(measures_of_fit);

% output table F test
disp(' ');
disp('F-test:');

F_test_table = round([ols.F_statistic; ols.p_valueF],6);

% convert to table
F_test_table = array2table(F_test_table);

% set row names
F_test_table.Properties.RowNames = {'F-statistic', 'F-stat p_value'};

% set column names
F_test_table.Properties.VariableNames = {' '};

% print result
disp(F_test_table);


% output table information criteria
disp(' ');
disp(' ');
disp('Criteria Model Selection:');

% array
info_criteria = round([ols.info_criteria.aic; ols.info_criteria.schwarz],4);

% convert to table
info_criteria = array2table(info_criteria);

% set row names
info_criteria.Properties.RowNames = {'AIC','Schwarz'};

% set column name
info_criteria.Properties.VariableNames = {' '};

% print results                                   
disp(info_criteria);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We hereby certify that
%We have written the program ourselves except for clearly marked pieces of code
%We have tested the program and it ran without crashing (if applicable)
%Matthias Lukosch (20-601-050)
%Chung Shun Man (20-621-587)
end