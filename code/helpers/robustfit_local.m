function varargout = robustfit_local(X,y,wfun,tune,const,priorw,dowarn)
%ROBUSTFIT Robust linear regression
%   B = ROBUSTFIT(X,Y) returns the vector B of regression coefficients,
%   obtained by performing robust regression to estimate the linear model
%   Y = Xb.  X is an n-by-p matrix of predictor variables, and Y is an
%   n-by-1 vector of observations.  The algorithm uses iteratively
%   reweighted least squares with the bisquare weighting function.  By
%   default, ROBUSTFIT adds a column of ones to X, corresponding to a
%   constant term in the first element of B.  Do not enter a column of ones
%   directly into the X matrix.
%
%   B = ROBUSTFIT(X,Y,'WFUN',TUNE) uses the weighting function 'WFUN' and
%   tuning constant TUNE.  'WFUN' can be any of 'andrews', 'bisquare',
%   'cauchy', 'fair', 'huber', 'logistic', 'talwar', or 'welsch'.
%   Alternatively 'WFUN' can be a function that takes a residual vector as
%   input and produces a weight vector as output.  The residuals are scaled
%   by the tuning constant and by an estimate of the error standard
%   deviation before the weight function is called.  'WFUN' can be
%   specified using @ (as in @myfun).  TUNE is a tuning constant that is
%   divided into the residual vector before computing the weights, and it
%   is required if 'WFUN' is specified as a function.
% 
%   B = ROBUSTFIT(X,Y,'WFUN',TUNE,'CONST') controls whether or not the
%   model will include a constant term.  'CONST' is 'on' (the default) to
%   include the constant term, or 'off' to omit it.
%
%   [B,STATS] = ROBUSTFIT(...) also returns a STATS structure
%   containing the following fields:
%       'ols_s'     sigma estimate (rmse) from least squares fit
%       'robust_s'  robust estimate of sigma
%       'mad_s'     MAD estimate of sigma; used for scaling
%                   residuals during the iterative fitting
%       's'         final estimate of sigma, the larger of robust_s
%                   and a weighted average of ols_s and robust_s
%       'se'        standard error of coefficient estimates
%       't'         ratio of b to stats.se
%       'p'         p-values for stats.t
%       'covb'      estimated covariance matrix for coefficient estimates
%       'coeffcorr' estimated correlation of coefficient estimates
%       'w'         vector of weights for robust fit
%       'h'         vector of leverage values for least squares fit
%       'dfe'       degrees of freedom for error
%       'R'         R factor in QR decomposition of X matrix
%
%   The ROBUSTFIT function estimates the variance-covariance matrix of the
%   coefficient estimates as V=inv(X'*X)*STATS.S^2.  The standard errors
%   and correlations are derived from V.
%
%   ROBUSTFIT treats NaNs in X or Y as missing values, and removes them.
%
%   Example:
%      x = (1:10)';
%      y = 10 - 2*x + randn(10,1); y(10) = 0;
%      bls = regress(y,[ones(10,1) x])
%      brob = robustfit(x,y)
%      scatter(x,y)
%      hold on
%      plot(x,brob(1)+brob(2)*x,'r-', x,bls(1)+bls(2)*x,'m:')
%
%   See also REGRESS, ROBUSTDEMO.

% References:
%   DuMouchel, W.H., and F.L. O'Brien (1989), "Integrating a robust
%     option into a multiple regression computing environment,"
%     Computer Science and Statistics:  Proceedings of the 21st
%     Symposium on the Interface, American Statistical Association.
%   Holland, P.W., and R.E. Welsch (1977), "Robust regression using
%     iteratively reweighted least-squares," Communications in
%     Statistics - Theory and Methods, v. A6, pp. 813-827.
%   Huber, P.J. (1981), Robust Statistics, New York: Wiley.
%   Street, J.O., R.J. Carroll, and D. Ruppert (1988), "A note on
%     computing robust regression estimates via iteratively
%     reweighted least squares," The American Statistician, v. 42,
%     pp. 152-154.

%   Copyright 1993-2011 The MathWorks, Inc.


if  nargin < 2      
    error(message('stats:robustfit:TooFewInputs'));      
end 

if (nargin<3 || isempty(wfun)), wfun = 'bisquare'; end
if nargin<4, tune = []; end
[wfun,tune] = statrobustwfun(wfun,tune);
if (nargin<5), const='on'; end
switch(const)
 case {'on' 1},  doconst = 1;
 case {'off' 0}, doconst = 0;
 otherwise,  error(message('stats:robustfit:BadConst'));
end
if nargin<6
    priorw = ones(size(y));
end
if nargin<7
    dowarn = true;
end

% Remove missing values and check dimensions
[anybad, wasnan, y, X, priorw] = removenan(y,X,priorw);
if (anybad==2)
   error(message('stats:robustfit:InputSizeMismatch'));
end

varargout=cell(1,max(1,nargout));
[varargout{:}] = statrobustfit_local(X,y,wfun,tune,wasnan,doconst,priorw,dowarn);



function [wfun,tune] = statrobustwfun(wfun,tune)
%STATROBUSTWFUN Get robust weighting function and tuning constant

% Copyright 2005-2007 The MathWorks, Inc.


% Convert name of weight function to a handle to a local function, and get
% the default value of the tuning parameter
t = [];
if ischar(wfun)
    switch(wfun)
        case 'andrews'
            wfun = @andrews;
            t = 1.339;
        case 'bisquare'
            wfun = @bisquare;
            t = 4.685;
        case 'cauchy'
            wfun = @cauchy;
            t= 2.385;
        case 'fair'
            wfun = @fair;
            t = 1.400;
        case 'huber'
            wfun = @huber;
            t = 1.345;
        case 'logistic'
            wfun = @logistic;
            t = 1.205;
        case 'ols'
            wfun = @ols;
            t = 1;
        case 'talwar'
            wfun = @talwar;
            t = 2.795;
        case 'welsch'
            wfun = @welsch;
            t = 2.985;
    end
end

% Use the default tuning parameter or check the supplied one
if isempty(tune)
    if isempty(t)
        tune = 1;
    else
        tune = t;
    end
elseif (tune<=0)
    m = message('stats:statrobustwfun:BadTuningConstant');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
     
% --------- weight functions

function w = andrews(r)
r = max(sqrt(eps(class(r))), abs(r));
w = (abs(r)<pi) .* sin(r) ./ r;

function w = bisquare(r)
w = (abs(r)<1) .* (1 - r.^2).^2;

function w = cauchy(r)
w = 1 ./ (1 + r.^2);

function w = fair(r)
w = 1 ./ (1 + abs(r));

function w = huber(r)
w = 1 ./ max(1, abs(r));

function w = logistic(r)
r = max(sqrt(eps(class(r))), abs(r));
w = tanh(r) ./ r;

function w = ols(r)
w = ones(size(r));

function w = talwar(r)
w = 1 * (abs(r)<1);

function w = welsch(r)
w = exp(-(r.^2));



function [badin,wasnan,varargout]=removenan(varargin)
%REMOVENAN Remove NaN values from inputs
%   [~,WASNAN,X1] = internal.stats.removenan(Y1) removes missing values
%   from Y1 and returns it as X1. WASNAN is a logical column vector
%   indicating where there were NaN values removed. Y1 is a column vector
%   or a matrix. The type of Y1 can be:
%
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [BAD,WASNAN,X1,x2,...] = internal.stats.removenan(Y1,Y2,...) accepts any
%  number of input variables Y1,Y2,Y3,..., removes missing values, and
%  returns them as X1, X2,... respectively. All variables should have the
%  same number of observations, and BAD is 0 if that is true. If they to
%  not have the same number, BAD returns the index of a variable whose
%  length does not match the others.
%
%  This utility is used by some Statistics Toolbox functions to handle
%  missing values.
%
%  See also INSERTNAN.

%   Copyright 1993-2012 The MathWorks, Inc.


badin = 0;
wasnan = 0;
n = -1;

% Find NaN, check length, and store outputs temporarily
varargout = cell(nargout,1);
for j=1:nargin
   y = varargin{j};
   if (size(y,1)==1) && (n~=1) 
       y =  y';
   end

   ny = size(y,1);
   if (n==-1)
      n = ny;
   elseif (n~=ny && ny~=0)
      if (badin==0), badin = j; end
   end
   
   varargout{j} = y;

   if (badin==0 && ny>0)
       wasnan = wasnan | any(isnan(y),2);
   end
end

if (badin>0), return; end

% Fix outputs
if (any(wasnan))
   t = ~wasnan;
   for j=1:nargin
      y = varargout{j};
      if (~isempty(y)), varargout{j} = y(t,:); end
   end
end
