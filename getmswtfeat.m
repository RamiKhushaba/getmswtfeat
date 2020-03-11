% GETMSWTFEAT Gets the Multiscale Wavelet Transform features, these
% include: Energy, Variance, Standard Deviation, and Waveform Length
% feat = getmswtfeat(x,winsize,wininc,SF)
% ------------------------------------------------------------------
% The signals in x are divided into multiple windows of size
% "winsize" and the windows are spaced "wininc" apart.
% Inputs
% ------
%    signals:	columns of signals
%    winsize:	window size (length of x)
%    wininc:	spacing of the windows (winsize)
%    SF:        sampling frequency (Not used in the current implementation, but I left you some options down there)
%
% Outputs
% -------
%    feat:     WT features organized as [Energy, Variance, Waveform Length, Entropy]

% Example
% -------
% feat = getmswtfeat(rand(1024,1),128,32,32)
% Assuming here rand(1024,1) (this can be any one dimensional signal,
% for example EEG or EMG) is a one dimensional signal sampled at 32
% for 32 seconds only. Utilizing a window size of 128 at 32 increments,
% features are extracted from the wavelet tree.
% I assumed 5 decomposition levels (J=5) below in the code.
% For a full tree at 5 levels you should get 6 features per window from each channel
% as we have decided to extract 5 types of features then we get 6 x 5 =30 features/channel
%
% =========================================================================
% Multiscale Wavelet Transform feature extraction code by Dr. Rami Khushaba
% Senior Research Fellow - ACFR
% University of Sydney
% Email: Rami.Khushaba@sydney.edu.au
% URL: www.rami-khushaba.com (Matlab Code Section)
% last modified 29/08/2012
% last modified 09/02/2013
% last modified 11/03/2020
% 
% P.S: the code is based on the template of the Myoelectric control
% developmenet toolbox by Adrian Chan: http://www.sce.carleton.ca/faculty/chan/index.php?page=matlab

function feature_out = getmswtfeat(signals,winsize,wininc,SF)

if nargin < 4
    if nargin < 3
        if nargin < 2
            error('A sliding window approach requires the window size (winsize) as input')
        end
        error('A sliding window approach requires the window increment (wininc) as input')
    end
    error('Please provide the sampling frequency of this signal')
end


%% The number of decomposition levels (manual or automatic?)
% Note I put SF above in the inputs because you can use SF to determine the
% best decompisition level J, however, for simplicity here I put it J=5;
decomOption = 1;

if decomOption==1
    J=5; % Number of decomposition levels set manually here
elseif decomOption==2
    J=wmaxlev(winsize,'Sym5'); % Number of decomposition levels set based on window size and wavelet family
else
    J=(log(SF/2)/log(2))-1; % Number of decomposition levels set based on sampling frequency (SF)
end

%% make sure you have some parameters pre-defined
% specify the number of samples
datasize = size(signals,1);
% based on the number of samples, winsize, and wininc, how many windows we
% will have? this is "numwin"
numwin = floor((datasize - winsize)/wininc)+1;
% how many signals (electrodes) are we processing
Nsignals = size(signals,2);
% how many features we plan to extract
NF = 5;
% predefine zeros matrix to allocate memory for output features
feature_out = zeros(numwin,(J+1)*NF*Nsignals);

for dims =1:Nsignals
    x=signals(:,dims);
    %% Chop the signal according to a sliding window approach
    % allocate memory
    feat = zeros(winsize,numwin);
    st = 1;
    en = winsize;
    for i = 1:numwin
        feat(1:winsize,i) = x(st:en,:)-mean(x(st:en,:));
        st = st + wininc;
        en = en + wininc;
    end
    
    %% Multisignal one-dimensional wavelet transform decomposition
    dec = mdwtdec('col',feat,J,'Sym8');
    % Proceed with Multisignal 1-D decomposition energy distribution
    
    if isequal(dec.dirDec,'c')
        dim = 1;
    end
    [cfs,longs] = wdec2cl(dec,'all');
    level = length(longs)-2;
    
    if dim==1
        cfs = cfs';
        longs = longs';
    end
    numOfSIGs             = size(cfs,1);
    num_CFS_TOT           = size(cfs,2);
    absCFS                = abs(cfs);
    absCFS0               = (cfs);
    cfs_POW2              = absCFS.^2;
    Energy                = sum(cfs_POW2,2);
    percentENER           = zeros(size(cfs_POW2));
    notZER                = (Energy>0);
    percentENER(notZER,:) = cfs_POW2(notZER,:);%./Energy(notZER,ones(1,num_CFS_TOT));
    
    %% or try this version below and tell us which  one is the best on your data
    % percentENER(notZER,:) = cfs_POW2(notZER,:);
    
    %% Pre-define and allocate memory
    tab_ENER    = zeros(numOfSIGs,level+1);
    tab_VAR     = zeros(numOfSIGs,level+1);
    tab_STD     = zeros(numOfSIGs,level+1);
    tab_WL      = zeros(numOfSIGs,level+1);
    tab_entropy = zeros(numOfSIGs,level+1);
    
    %% Feature extraction section
    st = 1;
    for k=1:level+1
        nbCFS                = longs(k);
        en                   = st+nbCFS-1;
        tab_ENER(:,k)        = sum(percentENER(:,st:en),2);% energy per waveform
        tab_VAR(:,k)         = var(percentENER(:,st:en),0,2); % variance of coefficients
        tab_STD(:,k)         = std(percentENER(:,st:en),[],2); % standard deviation of coefficients
        tab_WL(:,k)          = sum(abs(diff(percentENER(:,st:en)')))'; % waveform length
        percentENER(:,st:en) = percentENER(:,st:en)./repmat(sum(percentENER(:,st:en),2),1,size(percentENER(:,st:en),2));
        prob                 = percentENER(:,st:en)./repmat(sum(percentENER(:,st:en),2),1,longs(k)) + eps;
        tab_entropy(:,k)     = -sum(prob.*log(prob),2);%./size(percentENER(:,st:en),2);
        st = en + 1;
    end
    feature_out(:,(1:(NF*(J+1)))+(dims-1)*((J+1)*NF)) =log([tab_ENER tab_VAR tab_STD tab_WL tab_entropy]);
end