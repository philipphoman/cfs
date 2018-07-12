function out = cfs_phtdindiv(y, hybrid, lowerb, upperb, Xinit, fname);
% This function to calculate TD regressors and the hybrid pearce-hall
% model associability 
%      y:      string determining type of outcome
%      hybrid: string determining type of hybrid model 
%      lowerb: 1 x 4 vector of lower bounds
%      upperb: 1 x 4 vector of upper bounds
%      Xinit: starting values for X
%
%      out: summaries, 
%           HessSum, 
%           llss, 
%           TDs, 
%           alphaFits, 
%           VFits
%
% The update rule used here is: 
% alpha(n+1) = (1 - gamma)  * alpha(n) + gamma(abs(lambda(n) - V(n)))
% V(n+1)     = V(n) + kappa * alpha(n) * (lambda(n) - V(n))
%  
% However, Jian apparently used a softmax transformation of the term
% kappa * alpha 
%
%
%
% What the algorithm does: we want to find the set of free parameters
% (alpha, V, gamma, kappa) that maximize the likelihood function of our
% observed SCR data, or more formally:
% 
% P(SCR|M, theta_M), 
% 
% where M is comprised of (1) the learning model (i.e., the Pearce-Hall
% learning model) and (2) the observation model that tells us how the
% learning model's internal variables are reflected in observed data, 
% i.e., a linear regression model of the form 
%
% SCR_t = b0 + b1 * alpha_t + b2 * v_t + N(0, sigma)) 
% 
% and theta_M is the set of free parameters.
% 
% The task is to maximize the likelihood function for a given set of
% measured SCR. This can be done by searching for the full set of free
% parameters that optimize the likelihood function, which is what
% fmincon in the scripts in fact does. It iteratively tries different
% starting values for theta_M until the likelihood is maximized (which
% means that the negative log likelihood or the residual errors of (2)
% are minimized).
% 
% The optimal set of parameters is saved in the vector "XBest" in the
% script, the corresponding minimal value for the residual errors is
% saved as the scalar "tMin". The Hessian matrix for this solution is
% saved in "HessX" and is the second derivative of the likelihood
% function with respect to the parameters; it can be interpreted as the
% "charting" of a mountain, i.e it tells us how the hill of the data
% likelihood looks like, i.e., how much we can trust our estimate;
% larger values mean that the likelihood function is dropping off away
% from this estimate more rapidly, which corresponds to a more precise
% estimate. This is why we can use the inverse of the Hessian matrix
% (H^-1) as an estimator of the covariance of our estimate. While the
% diagonals in H^-1 are the variances of each parameter (that can be
% used to calculate standard errrors), the off-diagonals capture the
% covariance between the parameters, with large values indicating
% problems in the model or the data (so these should be useful for
% diagnostics).
% 
% Consider whether to apply restrictions on the free parameters. For
% fmincon, boundaries can be passed to the function for each
% parameter. Jian used [0 1] for alpha, v, and eta; and [0 inf] for
% kappa which intuitively makes sense. However, I read that boundaries
% should be used with caution as a true value may be lying outside the
% boundaries and may indicate problems with the model or a bad dataset
% where the observations aren't well explained by the model. Nathaniel
% suggests that one thing worthing trying is to use a prior instead of a
% hard constraint, i.e., using the population summary statistics (mean
% and variance) of the sample from which the subject was drawn as prior.
%-----------------------------------------------------------------------
% Original version by Jian Li
% Modified by PH, 1/2017
  
warning('off');
global cdat;
global tmpSCR;
global CSUS;
global US;
global hybrid;


if nargin < 5
  Xinit = nan(1, 4);
end

if nargin < 6
  fname       = '../data/ptsd.csv';
end

dat         = readDB(fname);
ntrials     = dat.ntrialsscr(1);
iterationN  = 20;
opt         = optimset('Display'     , 'off', ... 
                       'MaxFunEvals' , 10^5, ... 
                       'MaxIter'     , 10^5, ... 
                       'LargeScale'  , 'on', ... 
                       'Algorithm'   , 'active-set');
lowerBd     = lowerb;
higherBd     = upperb;

  

if (higherBd(4) == 1)
  suffix = '_allbounded';
else
  suffix = '';
end

%lowerBd     = [0 0 0 0];
%higherBd    = [1 1 1 1];
%higherBd    = [1 1 1 Inf];
%lowerBd     = [-Inf -Inf -Inf -Inf];
%higherBd    = [Inf Inf Inf Inf];
TDs         = [];
alphaFits   = [];
VFits       = [];
llss        = [];
BICs        = [];
SCRs        = [];
CSUSs       = [];
Orders      = [];
subj        = 1;

% discard subjects with fewer than ntrialsscr and with no dcm-data
if ~isempty(regexp(y, 'scr'))
  cdat   = dat(dat.ntrialsscr >= ntrials, :);
elseif ~isempty(regexp(y, 'dcm'))
  cdat   = dat(dat.ntrialsdcm >= ntrials, :);
end

% get subject list
idlist  = unique(cdat.id);

% cycle over subjects
%for k = 1:69:size(cdat, 1)
for k = 1:numel(idlist)
  idx = find(ismember(cdat.id, idlist{k}));
  subcdat = cdat(idx, :);
  subcdat = subcdat(1:ntrials, :);
%k = 1;
%while k < size(cdat, 1)
  fprintf(['\nProcessing subject #%i (%s) ', ...
           'starting with trial %i\n'], subj, subcdat.id{1},...
           subcdat.trial(1));

  %tmpSCR                 = log(subcdat.dcmadcm1 + 1);
  %tmpSCR                 = sqrt(subcdat.scr);
  %tmpSCR                 = subcdat.scrsqrtrc;
  tmpSCR                 = eval(['subcdat.' y]);
  CSUS                   = subcdat.stim;
  ord                    = subcdat.order;
  US                     = subcdat.us;
  X                      = randomX(Xinit);
  summaries{subj}        = [];
  HessSum{subj}          = [];

  [XBest, min, tmp2, ... 
   tmp3, tmp,tmp, HessX] = fmincon(@PHRegressors, ... 
                                   X, ... 
                                   [], ... 
                                   [], ... 
                                   [], ... 
                                   [], ... 
                                   lowerBd, ... 
                                   higherBd, ...
                                   [], ... 
                                   opt);

  summaries{subj}        = [summaries{subj}; [XBest min]];
  HessSum{subj}          = [HessSum{subj}; inv(-HessX)];
  %fprintf('Round 1 -->  The fitted parameters are:\n');
  fprintf('\tRound-%4d\t\t%12.3f\t%12.3f\t%12.3f\t%12.3f\n', 1, XBest);
  for i = 1:iterationN;
    X                       = randomX(Xinit);

    % save alternative Hessian as well
    [XBest2, tMin, tmp, ...
     tmp, tmp, tmp, HessX2] = fmincon(@PHRegressors, ... 
                                      X, ... 
                                      [], ... 
                                      [], ... 
                                      [], ... 
                                      [], ... 
                                      lowerBd, ... 
                                      higherBd, ... 
                                      [], ... 
                                      opt);
    
    summaries{subj}         = [summaries{subj}; [XBest2 tMin]];
    HessSum{subj}           = [HessSum{subj}; inv(-HessX2)];

    if (tMin < min & sum(XBest2 ~= XBest) > 0);
      XBest = XBest2;
      min   = tMin;

      % fprintf('\tRound-%4d\t\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n', i, XBest);
    end
  end

  % sort summaries according to minimal f(x) which is the last column
  summaries{subj}  = sortrows(summaries{subj}, size(lowerBd, 2)+1);
  
  % save first rows of summaries in extra matrix
  for s=1:numel(summaries), minsums(s,:) = summaries{s}(1, :); end

  % fit the free parameters to the observed data using optimized
  % starting values
  [resErrs, lls, TD, ... 
   alphaFit, VFit, bictmp] = PHRegressors(XBest);
  TDs              = [TDs; TD];
  alphaFits        = [alphaFits; alphaFit];
  VFits            = [VFits; VFit];

  % calculate log likelihoods for different values of gamma and kappa
  %lls              = [];
  llss             = [llss lls];
  
  % BIC = log(n) * k - 2 * log(L_hat)
  % k = q + 2 parameters, q is the number of slopes
  % Verified that this calculation is accurate
  BIC              = log(ntrials) * 2 - 2 * -lls;
  %[bictmp]
  BICs             = [BICs BIC];
  subj             = subj + 1;
  SCRs             = [SCRs tmpSCR];
  CSUSs            = [CSUSs CSUS];
  Orders           = [Orders ord];
end

% format and save output
%out              = struct('TDs'      , TDs,...
%                          'alphaFits', alphaFits,...
%                          'VFits'    , VFits,...
%                          'lls'      , lls,...
%                          'llss'     , llss,...
%                          'ids'      , cdat.id,...
%                          'summaries', summaries,...
%                          'HessSum'  , HessSum);
out              = [];
subids           = cdat.id;
%subidswide       = subids(1:69:end);
subidswide       = idlist;
save(['../data/cfs_phtdindiv_' hybrid '_' y suffix '.mat'], 'TDs', ...
     'alphaFits', 'VFits', 'lls', 'llss', 'idlist', 'summaries', ...
     'HessSum', 'minsums', 'BICs', 'SCRs', 'CSUSs', 'Orders')

function [resErrs, lls, TDs, alphaFits, VFits, bictmp] = PHRegressors(X);
  global cdat;
  global tmpSCR;
  global CSUS;
  global US;
  global hybrid;
  initCPAlpha = X(1); % initial associability for CS+ 
  initCMAlpha = X(1); % initial associability for CS-
  %initCPV     = X(2); % initial associative strength for CS+ 
  initCPV     = 0.5;
  %initCMV     = X(2); % initial associative strength for CS-
  initCMV     = 0.5;
  %eta         = X(3); % weight allocating factor
  eta         = 0;
  %kappa       = X(4); % scaling factor
  kappa       = 1;
  TDs         = [];
  alphaFits   = [];
  CSpTDs      = [];
  CSmTDs      = [];
  VFits       = [];
  bictmp      = [];
  CSpFirstInd = 0;
  CSmFirstInd = 0;
  initCVAlpha = X(1);
  initCV      = 0.5;
  CSUStmp     = CSUS;
  CSUStmp(2:end+1) = CSUS;
  for h = 1:length(CSUS)
    lambdabin(h) = strcmp(CSUS{h}, CSUStmp{h});
    lambda(h) = ~isempty(regexp(US{h}, 'US'));
  end
  lambdabin = abs(lambdabin - 1);
  lambdabin(1) = 0;
  for i = 1:size(CSUS, 1) 
    if strcmp(CSUS{i, 1}(1:6), 'CSplus')
      if ~CSpFirstInd;
        CPAlpha     = initCPAlpha;
        CPV         = initCPV;
        CSpFirstInd = 1;
        CVtmp       = initCV;
        CVAlpha     = initCPAlpha;
      else
        switch hybrid
          case 'rw'
            CPV         = CPV + CPAlpha * CSpTDs(end);

          case {'hybridv', 'hybridalpha', 'hybridalphav'}
            CPV         = CPV + kappa * CPAlpha * CSpTDs(end);
            %CPV         = CPV + (1/(1 + exp(-kappa * (CPAlpha -.5)))) * ...
            %             CSpTDs(end);

            CPAlpha     = (1 - eta) * CPAlpha + ...
                eta * abs(CSpTDs(end));
          
          case 'trialswitch'
            CVDelta     = lambdabin(i) - CVtmp;
            CVtmp       = CVtmp + CVAlpha * CVDelta;
            if CVtmp > 0.5
              %CV = CVtmp * round(CVtmp) + (1 - CVtmp) * (1 - round(CVtmp)); 
              s = abs(lambda(i-1)-1);
              CV = CVtmp * s + (1 - CVtmp) * (1 - s);
            else
              CV = CVtmp;
            end
            CPV         = CV;

        end
        

      end
      TDs       = [TDs;    strcmp(US{i, 1}(7:end), 'US') - CPV];
      CSpTDs    = [CSpTDs; strcmp(US{i, 1}(7:end), 'US') - CPV];
      alphaFits = [alphaFits; CPAlpha];
      VFits     = [VFits; CPV];
      
    elseif strcmp(CSUS{i, 1}(1:7), 'CSminus');
      if ~CSmFirstInd;
        CMAlpha     = initCMAlpha;
        CMV         = initCMV;
        CVtmp       = initCV;
        CVAlpha     = initCMAlpha;
        CSmFirstInd = 1;
      else
        switch hybrid
          case 'rw'
            CMV         = CMV + CMAlpha * CSmTDs(end);

          case {'hybridv', 'hybridalpha', 'hybridalphav'}
            CMV         = CMV + kappa * CMAlpha * CSmTDs(end);
            %CPV         = CPV + (1/(1 + exp(-kappa * (CPAlpha -.5)))) * ...
            %             CSpTDs(end);

            CMAlpha     = (1 - eta) * CMAlpha + ...
                eta * abs(CSmTDs(end));

          case 'trialswitch'
            CVDelta     = lambdabin(i) - CVtmp;
            CVtmp       = CVtmp + CVAlpha * CVDelta;
            if CVtmp > 0.5
              s = abs(lambda(i-1)-1);
              CV = CVtmp * s + (1 - CVtmp) * (1 - s);
            else
              CV = CVtmp;
            end
            CMV         = CV;
        end
        
        
        
      end
      TDs       = [TDs;    strcmp(US{i, 1}(8:end), 'US') - CMV];
      CSmTDs    = [CSmTDs; strcmp(US{i, 1}(8:end), 'US') - CMV];
      alphaFits = [alphaFits; CMAlpha];
      VFits     = [VFits; CMV];
    end
  end 
  USPInd           =  find(cellfun(@(x)strcmp(x,'CSplusUS'), ...
                                   CSUS(:,1)));
  USPInd           = [USPInd
                      find(cellfun(@(x)strcmp(x,'CSminusUS'), ...
                                   CSUS(:,1)))];

  USNInd           = setdiff([1:size(CSUS, 1)]', USPInd);
  alphaFitss       = alphaFits(USNInd);
  VFitss           = VFits(USNInd);
  rawSCRs          = tmpSCR(USNInd);
  
 
  switch hybrid
    case 'hybridalpha'
      [betas, resErrs, stats] = glmfit(alphaFitss, rawSCRs, ... 
                                'normal', 'link', 'identity');

      lls                     = -nansum(log(normpdf(rawSCRs, glmval(betas, ...
                                                             alphaFitss,...
                                                             'identity'), ...
                                                    stats.sfit)));

    case 'hybridalphav'
      [betas, resErrs, stats] = glmfit([alphaFitss VFitss], rawSCRs, ... 
                                'normal', 'link', 'identity');
      
      lls                     = -nansum(log(normpdf(rawSCRs, glmval(betas, ...
                                                             [alphaFitss ...
                                                              VFitss], ...
                                                             'identity'), ...
                                                    stats.sfit)));

    case {'rw', 'hybridv'}
      [betas, resErrs, stats] = glmfit(VFitss, rawSCRs, ... 
                                'normal', 'link', 'identity');
      
      
      %statstmp = fitglm(VFitss, rawSCRs);
      %bictmp = statstmp.ModelCriterion.BIC;

      lls                     = -nansum(log(normpdf(rawSCRs, glmval(betas, ...
                                                             VFitss, ...
                                                             'identity'), ...
                                                    stats.sfit)));
    case 'trialswitch'
      [betas, resErrs, stats] = glmfit(VFitss, rawSCRs, ... 
                                'normal', 'link', 'identity');

      

      %statstmp = fitglm(VFitss, rawSCRs);
      %bictmp = statstmp.ModelCriterion.BIC;
      
      lls                     = -nansum(log(normpdf(rawSCRs, glmval(betas, ...
                                                             VFitss, ...
                                                             'identity'), ...
                                                    stats.sfit)));
      
 
  end
  %lls              = -length(alphaFitss) * ...
  %                   (log(sqrt(2 * pi * resErrs/length(alphaFitss)))...
  %                    + 1/2);

return;


% read stimuli sequence from the database
function dat = readDB(fname);
  % fname                 = 'ptsd_phexport.xlsx';
  dat = readtable(fname);
return;

% Intitialize a starting point for all parameters
function X = randomX(Xinit);
  X = Xinit;
  X(isnan(Xinit)) = rand(1, sum(isnan(Xinit)));
  %X = [rand(1, 4)];
return;

% non-linear equality and inequality constraints
function [c, ceq] = confun(X);
  c   = [X(1) * X(4) - 1];
  ceq = [];  
return;
