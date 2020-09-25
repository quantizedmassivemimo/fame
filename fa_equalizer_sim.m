% =========================================================================
% -- Simulator for Finite-Alphabet Equalization
% -------------------------------------------------------------------------
% -- Simulator for the paper :
% -- Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein, and 
% -- Christoph Studer, "Finite-Alphabet MMSE Equalization for All-Digital
% -- Massive MU-MIMO mmWave Communication," IEEE Journal on Selected Areas
% -- in Communications (J-SAC), vol. 38, no. 9, September 2020
% -------------------------------------------------------------------------
% -- (c) 2020 Oscar Castañeda and Christoph Studer
% -- e-mail: oc66@cornell.edu and studer@cornell.edu
% =========================================================================

function fa_equalizer_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.runId = 0;       % simulation ID (used to reproduce results)    
    par.U = 16;          % number of single-antenna users
    par.B = 256;         % number of base-station antennas (B>>U)
    par.mod = '16QAM';   % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
    par.trials = 1e3;    % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = ... % list of signal-to-noise ratio [dB] values
        0:2:20;          % to be simulated
    par.equalizer = ...  % equalizing scheme(s) to be evaluated: 'LMMSE',
        ...              % 'FLMMSE','FAME_EXH_1b','FAME_SDR_1b','FAME_FBS'
        {'LMMSE','FLMMSE','FAME_FBS'};
    par.save = true;     % save results (true,false)
    par.plot = true;     % plot results (true,false)
    
    par.FAM.levels = 2^1; % number of levels to use in finite-alphabet
                          % matrices. A value of 2^B means that B bits are
                          % being completely used
            
    % *** FAME-FBS specific ***
    %
    % below we list the initial values used to train FAME-FBS using the 
    % 'deep unfolding' method described in the paper:
    %
    % Alexios Balatsoukas-Stimming, Oscar Castañeda, Sven Jacobsson, 
    % Giuseppe Durisi, and Christoph Studer, "Neural-Network Optimized
    % 1-Bit Precoding for Massive MU-MIMO," in Proceedings of the Inter-
    % national Workshop on Signal Processing Advances on Wireless Communi-
    % cations (SPAWC), July 2019
    %
    % BxU    | mod.  | levels | tau    | nu    | gamma | iters |
    % -------+-------+--------+--------+-------+-------+-------+
    % 8x2    | QPSK  | 2      | 2^(-4) | 1.25  | 1.1   | 30    |
    % -------+-------+--------+--------+-------+-------+-------+
    % 64x4   | 16QAM | 2      | 2^(-6) | 1.2   | 1.1   | 30    |
    % -------+-------+--------+--------+-------+-------+-------+
    % 256x16 | 16QAM | 2      | 2^(-8) | 1.15  | 1.1   | 5     |
    % -------+-------+--------+--------+-------+-------+-------+
    % 256x16 | 16QAM | 4      | 2^(-8) | 1.05  | 1.1   | 5     |
    % -------+-------+--------+--------+-------+-------+-------+
    % 256x16 | 16QAM | 8      | 2^(-8) | 1.05  | 1.1   | 5     |
    % -------+-------+--------+--------+-------+-------+-------+
    %
    % we provide the trained parameters for the configurations listed on
    % the previous table. To use these trained parameters, just set the
    % the number of FAME-FBS iterations to the number indicated by the
    % previous table, according to the system configuration being analyzed:
    
    par.FAME_FBS.iters = 5;
    
    % for other system configurations and FAME-FBS number of iterations,
    % please optimize the FAME-FBS parameters for best performance. For a
    % quick evaluation, you can start by keeping the FAME-FBS parameters
    % fixed for all iterations using the following:
    par.FAME_FBS.tau_fixed = 2^(-8);
    par.FAME_FBS.nu_fixed = 1.15;
    par.FAME_FBS.gamma_fixed = 1.1;
    % for the case were trained parameters are unavailable, you can also
    % select whether the low-resolution matrix X^H for FAME-FBS should be
    % initialized with the X^H from FL-MMSE or with the MRC solution H^H,
    % by setting the next variable to true of false, respectively:
    par.FAME_FBS.FLMMSEinit = false;
            
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1};   % only argument is par structure
    
end

%% Initialization

% set quantization parameters for channel matrix H and received vector y
id_str = [ 'B' num2str(par.B) '_U' num2str(par.U) '_Mod' par.mod ];
ids_quant = {'B8_U2_ModQPSK','B64_U4_Mod16QAM','B256_U16_Mod16QAM'};
if ismember(id_str,ids_quant)
    par.H_quant_bits = 8;
    par.y_quant_bits = 7;
    if strcmp(id_str,'B256_U16_Mod16QAM')
        par.H_quant_bits_int = 2;
        par.y_quant_bits_int = 5;
    else        
        par.H_quant_bits_int = 1;
        if strcmp(id_str,'B64_U4_Mod16QAM')
            par.y_quant_bits_int = 4;
        else % strcmp(id_str,'B8_U2_ModQPSK')
            par.y_quant_bits_int = 2;
        end
    end    
    par.H_y_quant = true;    
else
    warn_msg = [ 'Untested system configuration: Channel matrix H and ' ...
                 'received vector y will not be quantized. You can ' ...
                 'tune the quantization parameters if you wish.'];
    warning(warn_msg);
    par.H_y_quant = false;
end


% load trained parameters for FAME-FBS
if ismember('FAME_FBS',par.equalizer)
    par.FAME_FBS.tau = zeros(1,par.FAME_FBS.iters);
    par.FAME_FBS.nu = zeros(1,par.FAME_FBS.iters);
    par.FAME_FBS.gamma = zeros(1,par.FAME_FBS.iters);
    % - check if the file with parameters exists
    paramFileName = ['./params/famefbs_', ...
                     'B' num2str(par.B),'_U',num2str(par.U),'_',...
                     num2str(par.FAM.levels),'Levels_',...
                     num2str(par.FAME_FBS.iters), 'Iters.mat']; 
    if isfile(paramFileName)
        load(paramFileName,...
             'fbs_FLMMSEinit','fbs_tau','fbs_nu','fbs_gamma');
        par.FAME_FBS.FLMMSEinit = fbs_FLMMSEinit;
        par.FAME_FBS.tau(1,:) = fbs_tau;
        par.FAME_FBS.nu(1,:) = fbs_nu;
        par.FAME_FBS.gamma(1,:) = fbs_gamma;
    else
        warn_msg = [ 'No trained FAME-FBS parameters available for ' ...
                     'this system and FAME_FBS configuration. Using ' ...
                     'parameters in par.FAME_FBS. Note that parameter ' ...
                     'tuning is needed to obtain the best performance.' ];
        warning(warn_msg);
        par.FAME_FBS.tau(1,:) = ...
            par.FAME_FBS.tau_fixed * ones(1,par.FAME_FBS.iters);
        par.FAME_FBS.nu(1,:) = ...
            par.FAME_FBS.nu_fixed * ones(1,par.FAME_FBS.iters);
        par.FAME_FBS.gamma(1,:) = ...
            par.FAME_FBS.gamma_fixed * ones(1,par.FAME_FBS.iters);
    end
end

% use runId random seed (enables reproducibility)
rng(par.runId);

% simulation name (used for saving results)
par.simName = ['ERR_B',num2str(par.B),'_U',num2str(par.U), '_', ...
               par.mod, '_', num2str(par.trials),'Trials'];

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK'
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
end

% compute symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.bps = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% - initialize result arrays (equalizer x SNR)
% vector error rate
res.VER = zeros(length(par.equalizer),length(par.SNRdB_list));
% symbol error rate
res.SER = zeros(length(par.equalizer),length(par.SNRdB_list));
% bit error rate
res.BER = zeros(length(par.equalizer),length(par.SNRdB_list));
% denominator of error-vector magnitude (EVM)
res.VM = zeros(length(par.equalizer),length(par.SNRdB_list));
% mean-square error (note that this is across several channel realizations)
res.MSE = zeros(length(par.equalizer),length(par.SNRdB_list));
% SINDR
res.SINDR = zeros(length(par.equalizer),length(par.SNRdB_list));
% transmit power
res.TxPower = zeros(length(par.equalizer),length(par.SNRdB_list));
% receive power
res.RxPower = zeros(length(par.equalizer),length(par.SNRdB_list));
%noise power
res.NPower = zeros(length(par.equalizer),length(par.SNRdB_list));
% simulation equalization time
res.time = zeros(length(par.equalizer),length(par.SNRdB_list));

% compute noise variances to be considered
N0_list = par.U*par.Es*10.^(-par.SNRdB_list/10);

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.bps,par.trials);

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    x = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix and noise vector
    n = sqrt(0.5)*(randn(par.B,1)+1i*randn(par.B,1));
    H = sqrt(0.5)*(randn(par.B,par.U)+1i*randn(par.B,par.U));
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % set noise variance
        N0 = N0_list(k);
            
        % transmit data over noisy channel
        Hx = H*x;
        y = Hx + sqrt(N0)*n;
               
        % quantize entries of channel
        if par.H_y_quant
            HQ = QuantSatTrc(H,par.H_quant_bits_int, ...
                             par.H_quant_bits-par.H_quant_bits_int-1);
        else
            HQ = H;
        end
        
        % quantize received signals: models low-resolution ADCs at the
        % base-station
        if par.H_y_quant
            yQ = QuantSatTrc(y,par.y_quant_bits_int, ...
                             par.y_quant_bits-par.y_quant_bits_int-1);
        else
            yQ = y;
        end
        
        % algorithm loop
        for d=1:length(par.equalizer)                    
            
            % record time used by the equalizer
            starttime = toc;
            
            % generate equalizing matrix
            switch (par.equalizer{d})                
                case 'LMMSE'        % L-MMSE equalization
                    W_H = LMMSE(par, HQ, N0);
                case 'FLMMSE'       % FL-MMSE equalization
                    W_H = FLMMSE(par, HQ, N0);
                case 'FAME_EXH_1b'  % FAME via Exhaustive search
                    W_H = FAME_EXH_1b(par, HQ, N0);
                case 'FAME_SDR_1b'  % FAME via Semidefinite Relaxation
                    srng = rng;
                    [W_H] = FAME_SDR_1b(par,HQ,N0);
                    rng(srng);
                case 'FAME_FBS'     % FAME via FBS
                    W_H = FAME_FBS(par,HQ,N0);                
                otherwise
                    error('par.equalizer not specified')
            end
            
            % record equalization simulation time
            res.time(d,k) = res.time(d,k) + (toc-starttime);  
            
            %Apply equalization matrix
            xhat = W_H*yQ;
            
            % extract transmit and receive power
            res.TxPower(d,k) = res.TxPower(d,k) + ...
                               mean(sum(abs(x).^2))/par.U;
            res.RxPower(d,k) = res.RxPower(d,k) + ...
                               mean(sum(abs(Hx).^2));
            res.NPower(d,k) =  res.NPower(d,k) + ...
                               mean(sum(abs(sqrt(N0)*n).^2));
                                    
            % perform detection
            [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols)) ...
                         -ones(par.U,1)*par.symbols).^2,[],2);
            bithat = par.bits(idxhat,:);
            
            % -- compute error and complexity metrics
            err = (idx~=idxhat);
            res.VER(d,k) = res.VER(d,k) + any(err);
            res.SER(d,k) = res.SER(d,k) + sum(err)/par.U;
            res.BER(d,k) = res.BER(d,k) + ...
                           sum(sum(bits(:,:,t)~=bithat))/(par.U*par.bps);
            res.MSE(d,k) = res.MSE(d,k) + norm(xhat-x)^2;
            res.VM(d,k) = res.VM(d,k) + norm(x)^2;
            res.SINDR(d,k) = res.SINDR(d,k) + norm(x)^2/norm(xhat-x)^2;
            
        end % algorithm loop
        
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',...
            time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.EVM = sqrt(res.MSE./res.VM).*100;
res.MSE = res.MSE/par.trials;
res.SINDR = res.SINDR/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.NPower = res.NPower/par.trials;
res.time = res.time/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structures)

if par.save
    [~,~]=mkdir('./results');
    save(['./results/' par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plots)

if par.plot
    
    % - BER results
    marker_style = {'k-','b-','r-.','y-.','g-.','ms--','bv--'};
    figure()
    for d=1:length(par.equalizer)
        semilogy(par.SNRdB_list,res.BER(d,:), ...
                 marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    hold off
    grid on
    box on
    xlabel('average SNR per receive antenna [dB]','FontSize',12)
    ylabel('uncoded bit error-rate (BER)','FontSize',12);
    if length(par.SNRdB_list) > 1
        axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1]);
    end
    legend(par.equalizer, ...
           'FontSize',12,'location','southwest','Interpreter','None')
    set(gca,'FontSize',12);
    
    % - EVM results
    figure()
    for d=1:length(par.equalizer)
        plot(par.SNRdB_list,res.EVM(d,:),marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    % Mark minimum EVM (%) value required by 3GPP-5G NR in section 6.5.2.2
    minEVMmods = {'QPSK','16-QAM','64-QAM','256-QAM'};
    minEVM = [17.5, 12.5, 8, 3.5]; %{QPSK, 16QAM, 64QAM, 256QAM}
    for dd=1:length(minEVMmods)
        plot(par.SNRdB_list,minEVM(dd)*ones(size(par.SNRdB_list)), ...
             'r--','LineWidth',2);
        text(min(par.SNRdB_list), minEVM(dd)+1, minEVMmods{dd}, ...
             'Color','red');
    end
    hold off
    grid on
    box on
    xlabel('average SNR per receive antenna [dB]','FontSize',12)
    ylabel('error-vector magnitude (EVM) [%]','FontSize',12);
    if length(par.SNRdB_list) > 1
        axis([min(par.SNRdB_list) max(par.SNRdB_list) 0 40]);
    end
    legend(par.equalizer, ...
           'FontSize',12,'location','northeast','Interpreter','none')
    set(gca,'FontSize',12);
    
end

end

%% L-MMSE: Linear Minimum Mean-Square Error equalizer
function W_H = LMMSE(par,H,N0)

% compute equalization matrix with L-MMSE
W_H = (N0/par.Es*eye(par.U)+H'*H)\H';

end

%% FL-MMSE: Finite-alphabet Linear Minimum Mean-Square Error equalizer
function [V_H,X_H] = FLMMSE(par,H,N0)

% compute equalization matrix with L-MMSE
W_H = LMMSE(par,H,N0);

% quantize L-MMSE equalizer to obtain X^H
% - find largest real-valued entry of L-MMSE equalizer
W_R = real(W_H);
W_I = imag(W_H);
max_abs_W = max(abs([W_R W_I]),[],2);
max_abs_W = max_abs_W * ones(1,par.B);
% - quantizer on a per-row basis
% -- first, map to [-1,+1] range while saturating the extreme bins to
%    allow quantization using the rounding function
scale = 1-1/par.FAM.levels;
X_R = W_R./(scale*max_abs_W);
X_I = W_I./(scale*max_abs_W);
X_R(abs(X_R)>1) = sign(X_R(abs(X_R)>1));
X_I(abs(X_I)>1) = sign(X_I(abs(X_I)>1));
% -- map from [-1,+1] range to [0,1] to quantize with rounding function
X_R = round(0.5*(par.FAM.levels-1)*(X_R+1));
X_I = round(0.5*(par.FAM.levels-1)*(X_I+1));
% -- finish quantization by returning from [0,1] range to [-1,+1].
%    note that, following the paper notation, +1 corresponds to 
%    (1-1/par.FAM.levels)*w_max, but this scale does not need to be
%    recovered, as it will be handled by beta
X_R = (2/(par.FAM.levels-1))*X_R-1;
X_I = (2/(par.FAM.levels-1))*X_I-1;
X_H = X_R + 1i*X_I;

% compute scaling factors beta to complete the finite-alphabet matrix V^H
V_H = zeros(size(X_H));
for uu=1:par.U
    hu = H(:,uu);
    xu_h = X_H(uu,:);
    beta = xu_h*hu/(norm(H'*xu_h')^2+N0/par.Es*norm(xu_h)^2);
    vu_h = beta'*xu_h;
    V_H(uu,:) = vu_h;
end

end


%% FAME-EXH: Finite-Alphabet Minimum mean-square Error equalizer solved
%            via exhaustive search, for 1-bit alphabets
function V_H = FAME_EXH_1b(par,H,N0)

if par.B>8
    err_msg = ['You are running an exhaustive search for a ' ...
               'base-station with more than 8 antennas. ' ...
               'This will take a LONG time. ' ...
               'Do you really want to proceed?'];
    error(err_msg);
end

if par.FAM.levels~=2
    err_msg = ['The provided FAME-EXH function only works for 1-bit ' ...
               'alphabets. Please set par.FAM.levels=2 or remove ' ...
               'FAME_EXH_1b from par.equalizer.'];
    error(err_msg);
end

rho = N0/par.Es;

% note that this code is for 1-bit alphabets only
symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
symbols_comb = de2bi((1:4^par.B)-1,par.B,4)'+1;
X = symbols(symbols_comb);

% compute all possible quantities for FAME cost function
numerator = sum(abs(H'*X).^2,1) + rho*sum(abs(X).^2,1);
denominators = abs(H'*X).^2;

V_H = zeros(par.U,par.B);
for uu=1:par.U
    % find vector x that minimizes the cost function
    cost = numerator./denominators(uu,:);
    [~,idx] = min(cost);
    x = X(:,idx);
    % compute corresponding scaling factor beta
    beta = (x'*H(:,uu))/numerator(idx);
    % compute finite-alphabet matrix row
    vu = beta*X(:,idx);
    V_H(uu,:) = vu';
end

end

%% FAME-SDR: Finite-Alphabet Minimum mean-square Error equalizer solved
%            using semidefinite-relaxation (SDR), for 1-bit alphabets
function V_H = FAME_SDR_1b(par,H,N0)

if par.FAM.levels~=2
    err_msg = ['The provided FAME-EXH function only works for 1-bit ' ...
               'alphabets. Please set par.FAM.levels=2 or remove ' ...
               'FAME_EXH_1b from par.equalizer.'];
    error(err_msg);
end

% convert to real-valued channel
HR = [real(H), -imag(H); imag(H), real(H)];

% set up SDP
rho = N0/par.Es;
T = HR*HR'+rho*eye(par.B*2);

V_H = zeros(par.U,par.B);

for uu=1:par.U
    
    huR =  HR(:,uu);

    % solve SDP
    cvx_begin quiet
        cvx_precision default
        variable X(2*par.B,2*par.B) symmetric
        minimize trace(T*X)
        subject to
            for jj = 2:2*par.B
                X(jj,jj) == X(1,1);
            end
            huR'*X*huR == 1;
            X == semidefinite(2*par.B);
    cvx_end

    % rank-one approximation
    % - eigenvalue decomposition
    [eigV, eigD] = eig(X);
    % - find maximum eigenvalue
    [~, idxmax] = max(diag(eigD));
    % - approximate solution
    xR = sqrt(eigD(idxmax,idxmax))*eigV(:,idxmax);
    xR = sign(xR);
    % - recover binary, complex-valued solution
    x = xR(1:par.B,1)+1i*xR(par.B+1:2*par.B,1);
    
    % compute corresponding scaling factor beta
    hu = H(:,uu);
    beta = (x'*hu)/(norm(H'*x)^2+rho*norm(x)^2);
    % compute finite-alphabet matrix row
    vu = beta*x;
    V_H(uu,:) = vu';
    
end

end

%% FAME-FBS: Finite-Alphabet Minimum mean-square Error equalizer solved
%            via Forward-Backward Splitting
function V_H = FAME_FBS(par,H,N0)

rho = N0/par.Es;

% choose initializer for low-resolution matrix X^H
if par.FAME_FBS.FLMMSEinit
    [~,X_H] = FLMMSE(par,H,N0);
    X = X_H';
else
    X = H;
end

V_H = zeros(par.U,par.B);

for uu=1:par.U
    
    % use FBS to find low-resolution matrix X^H
    eu_mat = zeros(par.U);
    eu_mat(uu,uu) = 1;    
    x = X(:,uu);
    for i=1:par.FAME_FBS.iters
        % - gradient descent, equation (32)
        Hhx = H'*x;
        diffIU = Hhx - par.FAME_FBS.gamma(i)*eu_mat*Hhx;
        x = x-par.FAME_FBS.tau(i)*H*diffIU;
        % - projection using proximal operator, equation (33) and (31)
        x = par.FAME_FBS.nu(i)*x;        
        x = min(max(real(x),-1),1) + 1i*min(max(imag(x),-1),1);
    end
    
    % quantize final solution of FBS
    x_R = real(x);
    x_I = imag(x);
    % - map from [-1,+1] range to [0,1] to quantize with rounding function    
    x_R = round(0.5*(par.FAM.levels-1)*(x_R+1));
    x_I = round(0.5*(par.FAM.levels-1)*(x_I+1));
    % - finish quantization by returning from [0,1] range to [-1,+1]
    x_R = (2/(par.FAM.levels-1))*x_R-1;
    x_I = (2/(par.FAM.levels-1))*x_I-1;
    x = x_R + 1i*x_I;
    
    % compute scaling factor beta
    hu = H(:,uu);    
    beta = (x'*hu)/(norm(H'*x)^2+rho*norm(x)^2);
    % compute finite-alphabet matrix row
    vu = beta*x;
    V_H(uu,:) = vu';
    
end

end

%% QuantSatTrc
% Quantize with saturation and truncation, for signed numbers (so the total
% number of bits needed is int_bits+frac_bits+1 (sign_bit)
% saturation is symmetrical (saturated at same value for both positive and
% negative numbers)
function XQ = QuantSatTrc(X,int_bits,frac_bits)

max_abs = (2^(int_bits+frac_bits)-1)/2^frac_bits;

% Quantize real part
X_R = real(X);
X_R(abs(X_R)>max_abs) = max_abs*sign(X_R(abs(X_R)>max_abs)); % saturate
XQ_R = floor(X_R*2^frac_bits)/2^frac_bits; % truncate

% Quantize imaginary part
X_I = imag(X);
X_I(abs(X_I)>max_abs) = max_abs*sign(X_I(abs(X_I)>max_abs)); % saturate
XQ_I = floor(X_I*2^frac_bits)/2^frac_bits; % truncate

XQ = XQ_R + 1j*XQ_I;

end