function [rmse, tgtposEst, tgtposTrue, anchorpos] = runLocalizationScenario(params)
%RUNLOCALIZATIONSCENARIO  Single TOA/TDOA localization scenario.
%
%   params struct fields:
%     numAnchors      (int)    number of anchors  (>=3)
%     numTargets      (int)    number of targets  (>=1)
%     measurementType (char)   'TOA' or 'TDOA'
%     spectrumMethod  (char)   'FFT' or 'MUSIC'
%     bw              (Hz)     bandwidth
%     Pt              (W)      transmit peak power
%     NF              (dB)     receiver noise figure
%     delayoffset     (s)      synchronization offset
%     N               (int)    OFDM subcarriers
%     M               (int)    OFDM symbols
%     seed            (int)    rng seed (controls anchor/target geometry & noise)
%
%   Returns RMSE (m) of estimated target positions, plus the raw
%   estimated/true target and anchor positions.

    rng(params.seed);

    fc = 38e9;
    cLight = physconst('LightSpeed');
    sampleRate = params.bw;
    Gtx = 20;
    Grx = 20;

    numAnchors = params.numAnchors;
    numTargets = params.numTargets;
    N = params.N;
    M = params.M;

    % Random anchor positions
    anchorpos = [60*rand(1, numAnchors) - 30; ...
                 60*rand(1, numAnchors) - 30; ...
                 40*rand(1, numAnchors) - 20];

    % Random target positions within anchor bounding box
    xRange = [min(anchorpos(1,:)), max(anchorpos(1,:))];
    yRange = [min(anchorpos(2,:)), max(anchorpos(2,:))];
    zRange = [min(anchorpos(3,:)), max(anchorpos(3,:))];
    tgtposTrue = [ ...
        xRange(1) + diff(xRange) * rand(1, numTargets); ...
        yRange(1) + diff(yRange) * rand(1, numTargets); ...
        zRange(1) + diff(zRange) * rand(1, numTargets)];

    % OFDM waveform parameters
    freqSpacing = params.bw / N;
    maxDelay    = 200e-9;
    rmax        = maxDelay * cLight;
    tcp         = range2time(rmax);
    Ncp         = ceil(sampleRate * tcp);
    tWave       = (1/freqSpacing) + Ncp/sampleRate;
    Ns          = N + Ncp;

    % Transceiver chain
    antenna     = phased.IsotropicAntennaElement('BackBaffled', false);
    transmitter = phased.Transmitter('Gain', Gtx, 'PeakPower', params.Pt);
    radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc);
    collector   = phased.Collector('Sensor', antenna, 'OperatingFrequency', fc);
    receiver    = phased.Receiver('AddInputNoise', true, 'Gain', Grx, ...
                  'NoiseFigure', params.NF, 'SampleRate', sampleRate);
    channel     = phased.FreeSpace('PropagationSpeed', cLight, ...
                  'OperatingFrequency', fc, 'SampleRate', sampleRate, ...
                  'TwoWayPropagation', false);

    tgtposEst = zeros(3, numTargets);

    for idxTgt = 1:numTargets
        tgtpos = tgtposTrue(:, idxTgt);
        tgtplatform    = phased.Platform('InitialPosition', tgtpos, 'Velocity', [0;0;0]);
        anchorplatform = phased.Platform('InitialPosition', anchorpos, 'Velocity', zeros(3, numAnchors));

        X = cell(1, numAnchors);
        for idxAnchor = 1:numAnchors
            bpskSymbol = randi([0,1], [N M]) * 2 - 1;
            sigmod = ifft(ifftshift(bpskSymbol, 1), [], 1);
            sig = sigmod([end-Ncp+(1:Ncp), 1:end], :);
            sig = sig / max(abs(sig), [], 'all');

            x = complex(zeros(size(sig)));
            for m = 1:M
                [tx_pos, tx_vel] = anchorplatform(tWave);
                [rx_pos, rx_vel] = tgtplatform(tWave);
                [~, txang] = rangeangle(rx_pos, tx_pos(:, idxAnchor));
                txsig = transmitter(sig);
                radtxsig = radiator(txsig(:, m), txang);
                chansig = channel(radtxsig, tx_pos(:, idxAnchor), rx_pos, ...
                    tx_vel(:, idxAnchor), rx_vel);
                [~, rxang] = rangeangle(tx_pos(:, idxAnchor), rx_pos);
                rxsig = collector(chansig, rxang);
                x(:, m) = receiver(rxsig);
            end

            xvec = reshape(x, Ns*M, 1);
            xdemod = ofdmdemod(xvec, N, Ncp, Ncp);
            X{idxAnchor} = xdemod ./ bpskSymbol;

            reset(anchorplatform);
            reset(tgtplatform);
        end

        if strcmp(params.spectrumMethod, 'FFT')
            estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                'Measurement', params.measurementType, 'SpectrumMethod', 'FFT', ...
                'VarianceOutputPort', true, 'DelayOffsetInputPort', true);
        else
            estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                'Measurement', params.measurementType, 'SpectrumMethod', 'MUSIC', ...
                'VarianceOutputPort', true, 'DelayOffsetInputPort', true, ...
                'ForwardBackwardAveraging', true, 'SpatialSmoothing', ceil(N/2));
        end

        [Y, estVar] = estimator(X, freqSpacing, params.delayoffset);
        if strcmp(params.measurementType, 'TOA')
            tgtposEst(:, idxTgt) = toaposest(Y, estVar, anchorpos);
        else
            tgtposEst(:, idxTgt) = tdoaposest(Y, estVar, anchorpos);
        end
    end

    errors = sqrt(sum((tgtposEst - tgtposTrue).^2, 1));
    rmse = sqrt(mean(errors.^2));
end
