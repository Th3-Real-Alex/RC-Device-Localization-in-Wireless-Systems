
%Scenario Configuration

rng('default');

% RF parameters
fc = 38e9;                                       % Carrier frequency (Hz)
c = physconst('LightSpeed');                     % Speed of propagation (m/s)
bw = 100e6;                                      % Bandwidth (Hz)
sampleRate = bw;                                 % Sample rate (Hz)

% Tx and Rx parameters
Pt = 0.1;                                        % Peak power (W) 
Gtx = 20;                                        % Tx antenna gain (dB)
Grx = 20;                                        % Rx antenna gain (dB) 
NF = 2.9;                                        % Noise figure (dB) 

% Create a transceiver components to simulate waveform propagation
antenna = phased.IsotropicAntennaElement('BackBaffled',false);
transmitter = phased.Transmitter('Gain',Gtx,'PeakPower',Pt);
radiator = phased.Radiator('Sensor',antenna,'OperatingFrequency',fc);
collector = phased.Collector('Sensor',antenna,'OperatingFrequency',fc);
receiver = phased.Receiver('AddInputNoise',true,'Gain',Grx, ...
         'NoiseFigure',NF,'SampleRate',sampleRate);

% Create a one-way free-space propagation channel
channel = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc, ...
    'SampleRate',sampleRate,'TwoWayPropagation',false);

% Device platform
tgtpos = [8; 4; 0];                              % Device position (m)
tgtvel = [1; -1; 0.5];                           % Device velocity (m/s)
tgtplatform = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel); 

% Anchor platform
anchorpos = [0 30 10 20 15; ...
    10 -5 -20 15 10; ...
    4 -10 10 5 -20];                             % Anchor positions (m)
numAnchor = size(anchorpos,2);                   % Number of anchors 
anchorvel = zeros(3,numAnchor);                  % Anchor velocities (m/s)
anchorplatform = phased.Platform('InitialPosition',anchorpos,'Velocity',anchorvel);







%Waveform Configuration
% Configure OFDM waveform
N = 1024;                                        % Number of subbands (OFDM data carriers)
M = 8;                                           % Number of channel samples (OFDM symbols)
freqSpacing = bw/N;                              % Frequency spacing (Hz)
tsym = 1/freqSpacing;                            % OFDM symbol duration  (s)
maxDelay = 200e-9;                               % Maximum delay (s)
rmax = maxDelay*c;                               % Maximum range of interest (m)
tcp = range2time(rmax);                          % Duration of the CP (s)
Ncp = ceil(sampleRate*tcp);                      % Length of the CP in subcarriers
tcp = Ncp/sampleRate;                            % Adjusted duration of the CP (s)
tWave = tsym + tcp;                              % OFDM symbol duration with CP (s)
Ns = N + Ncp;                                    % Number of subcarriers in one OFDM symbol
%Channel Estimation
% Estimated channel
X = cell(1,numAnchor);

% Flag to indicate if Communications Toolbox is licensed
useCommunicationsToolbox  = false;

% OFDM waveform propagation and channel estimation
for idxAnchor = 1:numAnchor
    % BPSK symbol on each data carrier
    bpskSymbol = randi([0,1],[N M])*2-1;             
    
    % Modulate OFDM waveform with BPSK modulation on each data carrier
    if useCommunicationsToolbox
        % OFDM modulated signal with cyclic prefix
        sigmod = ofdmmod(bpskSymbol,N,Ncp); %#ok<UNRCH>
        % Reshape OFDM modulated signal
        sig = reshape(sigmod,Ns,M);                  
    else
        % OFDM symbols
        sigmod = ifft(ifftshift(bpskSymbol,1),[],1);
        % OFDM symbols with cyclic prefix
        sig = sigmod([end-Ncp+(1:Ncp),1:end],:);     
    end

    % Power-normalized OFDM signal
    sig = sig/max(abs(sig),[],'all');

    % Initialize estimated channel
    x = complex(zeros(size(sig)));

    % Transceiver chain
    for m = 1:M
        % Update radar and device positions
        [tx_pos,tx_vel] = anchorplatform(tWave);
        [rx_pos,rx_vel] = tgtplatform(tWave);

        % Calculate the transmit angle
        [~,txang] = rangeangle(rx_pos,tx_pos(:,idxAnchor)); 

        % Transmitted signal
        txsig = transmitter(sig);

        % Radiate signal towards the receiver 
        radtxsig = radiator(txsig(:,m),txang); 

        % Propagate the signal
        chansig = channel(radtxsig,tx_pos(:,idxAnchor),rx_pos, ...
            tx_vel(:,idxAnchor),rx_vel);

        % Calculate the receive angle
        [~,rxang] = rangeangle(tx_pos(:,idxAnchor),rx_pos);

        % Collect signal at the receive antenna
        rxsig = collector(chansig,rxang);

        % Receive signal at the receiver
        x(:,m) = receiver(rxsig);
    end

    if useCommunicationsToolbox
        % Remove the cyclic prefix 
        xrmvcp = x((Ncp+1):(Ncp+N),:); %#ok<UNRCH>

        % Demodulate the received OFDM signal
        xdemod = fftshift(fft(xrmvcp,[],1),1);
    else
        xvec = reshape(x,Ns*M,1);

        % Remove the cyclic prefix and demodulate the received OFDM signal
        xdemod = ofdmdemod(xvec,N,Ncp,Ncp);
    end

    % Remove BPSK symbol on each data carrier
    X{idxAnchor} = xdemod./bpskSymbol;

    % Reset platforms for the next anchor
    reset(anchorplatform);
    reset(tgtplatform);
end






%TOA Estimation and Localization
% Spectrum analysis method
spectrumMethod = "FFT";
 
% Configure TOA estimator
if strcmp(spectrumMethod,'FFT') 
    toaEstimator = phased.TOAEstimator('PropagationSpeed',c, ...
        'Measurement','TOA','SpectrumMethod',spectrumMethod, ... 
        'VarianceOutputPort',true,'DelayOffsetInputPort',true); 
else % 'MUSIC'
    toaEstimator = phased.TOAEstimator('PropagationSpeed',c, ...
        'Measurement','TOA','SpectrumMethod',spectrumMethod, ... 
        'VarianceOutputPort',true,'DelayOffsetInputPort',true, ...
        'ForwardBackwardAveraging',true,'SpatialSmoothing',ceil(N/2)); %#ok<UNRCH>
end
% Synchronization error in seconds
delayoffset = 0;

% TOA estimation
[Y,var] = toaEstimator(X,freqSpacing,delayoffset);
% Plot TOA spectrum for anchor 1
figure
plotTOASpectrum(toaEstimator,freqSpacing,'AnchorIndex',1,'MaxDelay',maxDelay);
% Obtain TOA position estimate
tgtposest = toaposest(Y,var,anchorpos);
% View TOA position estimate
helperPlotTOADevicePositions(c,Y,tgtposest,anchorpos,tgtpos);


%TDOA Estimation and Localization
% Spectrum analysis method
spectrumMethod = "MUSIC";
 
% Configure TDOA estimator
if strcmp(spectrumMethod,'FFT') 
    tdoaEstimator = phased.TOAEstimator('PropagationSpeed',c, ...
        'Measurement','TDOA','SpectrumMethod',spectrumMethod, ... 
        'VarianceOutputPort',true,'DelayOffsetInputPort',true); %#ok<UNRCH>
else % 'MUSIC'
    tdoaEstimator = phased.TOAEstimator('PropagationSpeed',c, ...
        'Measurement','TDOA','SpectrumMethod',spectrumMethod, ... 
        'VarianceOutputPort',true,'DelayOffsetInputPort',true, ...
        'ForwardBackwardAveraging',true,'SpatialSmoothing',ceil(N/2)); 
end
% Common synchronization error in seconds
delayoffset = 100e-9;

% TDOA estimation
[Y,var] = tdoaEstimator(X,freqSpacing,delayoffset);
% Obtain TDOA position estimate
tgtposest = tdoaposest(Y,var,anchorpos);
% View TDOA position estimate
helperPlotTDOADevicePositions(c,Y,tgtposest,anchorpos,tgtpos)
% RMSE of the TDOA position estimate
RMSE = rmse(tgtposest,tgtpos);
disp(['RMS Localization error = ', num2str(RMSE), ' meters.'])