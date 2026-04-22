function positionAccuracy()
    pageContainer = uifigure('Name', 'Device Localization in Wireless Systems - Position Accuracy', ...
        'NumberTitle', 'off', ...
        'Resize', 'off', ...
        'ToolBar', 'none', ...
        'Position', [200 200 1280 720] ...
    );

    % ---- Constants ----
    panelWidth = 310;
    leftPad = 15;
    ctrlWidth = 280;
    ctrlHeight = 22;
    labelHeight = 16;
    rowSpacing = 34;

    % ---- Left Panel (controls) ----
    controlPanel = uipanel(pageContainer, ...
        'Title', 'Simulation Parameters', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'Position', [10 10 panelWidth 700] ...
    );

    % Build controls from top to bottom
    yPos = 640;

    % --- Measurement Type ---
    yPos = yPos - 5;
    uilabel(controlPanel, 'Text', 'Measurement Type', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddMeasurement = uidropdown(controlPanel, ...
        'Items', {'TOA', 'TDOA'}, ...
        'Value', 'TOA', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Spectrum Method ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Spectrum Method', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSpectrum = uidropdown(controlPanel, ...
        'Items', {'FFT', 'MUSIC'}, ...
        'Value', 'FFT', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Number of Anchors ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Number of Anchors (3-10)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnAnchors = uispinner(controlPanel, ...
        'Value', 5, 'Limits', [3 10], 'Step', 1, ...
        'RoundFractionalValues', 'on', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Number of Targets ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Number of Targets (1-5)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnTargets = uispinner(controlPanel, ...
        'Value', 1, 'Limits', [1 5], 'Step', 1, ...
        'RoundFractionalValues', 'on', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Bandwidth (MHz) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Bandwidth (MHz, 50-400)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnBandwidth = uispinner(controlPanel, ...
        'Value', 100, 'Limits', [50 400], 'Step', 50, ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Transmit Power (W) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Transmit Power (W, 0.01-1)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnPower = uispinner(controlPanel, ...
        'Value', 0.1, 'Limits', [0.01 1], 'Step', 0.05, ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Noise Figure (dB) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Noise Figure (dB, 1-10)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnNoise = uispinner(controlPanel, ...
        'Value', 2.9, 'Limits', [1 10], 'Step', 0.5, ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Delay Offset (ns) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Delay Offset (ns, 0-500)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnDelay = uispinner(controlPanel, ...
        'Value', 0, 'Limits', [0 500], 'Step', 10, ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Num Subcarriers (N) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Subcarriers (N)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSubcarriers = uidropdown(controlPanel, ...
        'Items', {'256', '512', '1024', '2048'}, ...
        'Value', '1024', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Num OFDM Symbols (M) ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'OFDM Symbols (M)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], ...
        'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSymbols = uidropdown(controlPanel, ...
        'Items', {'4', '8', '16'}, ...
        'Value', '8', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    % --- Run Button ---
    yPos = yPos - rowSpacing - 5;
    btnRun = uibutton(controlPanel, ...
        'Text', 'Run Simulation', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.93 0.53 0.24], ...
        'FontColor', 'white', ...
        'Position', [leftPad yPos ctrlWidth 35], ...
        'ButtonPushedFcn', @(~,~) runSimulation());

    % --- RMSE Display ---
    yPos = yPos - 35;
    lblRMSE = uilabel(controlPanel, ...
        'Text', 'RMS Error: —', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'FontColor', [0.2 0.2 0.8], ...
        'Position', [leftPad yPos ctrlWidth 30]);

    % ---- Right Panel (axes) ----
    ax = uiaxes(pageContainer, ...
        'Position', [panelWidth + 30, 30, 1280 - panelWidth - 50, 660]);
    title(ax, 'Device Localization Results');
    xlabel(ax, 'X (meters)');
    ylabel(ax, 'Y (meters)');
    grid(ax, 'on');
    hold(ax, 'on');

    % ---- Status Label (above axes) ----
    lblStatus = uilabel(pageContainer, ...
        'Text', 'Ready. Configure parameters and click Run Simulation.', ...
        'FontSize', 12, ...
        'FontColor', [0.3 0.3 0.3], ...
        'Position', [panelWidth + 30, 690, 600, 25]);

    % ======================================================================
    %  SIMULATION ENGINE (nested function — has access to all UI handles)
    % ======================================================================
    function runSimulation()
        % --- Disable Run button & show busy ---
        btnRun.Enable = 'off';
        btnRun.Text = 'Running...';
        lblStatus.Text = 'Simulation in progress...';
        lblRMSE.Text = 'RMSE: —';
        drawnow;

        try
            % --- Read parameters from UI ---
            measurementType = ddMeasurement.Value;
            spectrumMethod  = ddSpectrum.Value;
            numAnchors      = spnAnchors.Value;
            numTargets      = spnTargets.Value;
            bw              = spnBandwidth.Value * 1e6;         % Hz
            Pt              = spnPower.Value;                   % W
            NF              = spnNoise.Value;                   % dB
            delayoffset     = spnDelay.Value * 1e-9;            % s
            N               = str2double(ddSubcarriers.Value);
            M               = str2double(ddSymbols.Value);

            % --- Fixed RF parameters ---
            fc = 38e9;                                          % Carrier frequency (Hz)
            cLight = physconst('LightSpeed');
            sampleRate = bw;
            Gtx = 20;                                          % Tx antenna gain (dB)
            Grx = 20;                                          % Rx antenna gain (dB)

            % --- Generate random anchor positions ---
            rng('default');
            anchorpos = [60*rand(1, numAnchors) - 30; ...       % X: [-30, 30]
                         60*rand(1, numAnchors) - 30; ...       % Y: [-30, 30]
                         40*rand(1, numAnchors) - 20];          % Z: [-20, 20]

            % --- Generate random target positions (within anchor bounds) ---
            xRange = [min(anchorpos(1,:)), max(anchorpos(1,:))];
            yRange = [min(anchorpos(2,:)), max(anchorpos(2,:))];
            zRange = [min(anchorpos(3,:)), max(anchorpos(3,:))];
            tgtposAll = [ ...
                xRange(1) + diff(xRange) * rand(1, numTargets); ...
                yRange(1) + diff(yRange) * rand(1, numTargets); ...
                zRange(1) + diff(zRange) * rand(1, numTargets)];

            % --- OFDM waveform parameters ---
            freqSpacing = bw / N;
            maxDelay    = 200e-9;
            rmax        = maxDelay * cLight;
            tcp         = range2time(rmax);
            Ncp         = ceil(sampleRate * tcp);
            tWave       = (1/freqSpacing) + Ncp/sampleRate;
            Ns          = N + Ncp;

            % --- Create transceiver components ---
            antenna     = phased.IsotropicAntennaElement('BackBaffled', false);
            transmitter = phased.Transmitter('Gain', Gtx, 'PeakPower', Pt);
            radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc);
            collector   = phased.Collector('Sensor', antenna, 'OperatingFrequency', fc);
            receiver    = phased.Receiver('AddInputNoise', true, 'Gain', Grx, ...
                          'NoiseFigure', NF, 'SampleRate', sampleRate);
            channel     = phased.FreeSpace('PropagationSpeed', cLight, ...
                          'OperatingFrequency', fc, 'SampleRate', sampleRate, ...
                          'TwoWayPropagation', false);

            % --- Per-target simulation ---
            tgtposEstAll = zeros(3, numTargets);

            for idxTgt = 1:numTargets
                tgtpos = tgtposAll(:, idxTgt);
                tgtvel = [0; 0; 0];

                % Platforms
                tgtplatform    = phased.Platform('InitialPosition', tgtpos, 'Velocity', tgtvel);
                anchorvel      = zeros(3, numAnchors);
                anchorplatform = phased.Platform('InitialPosition', anchorpos, 'Velocity', anchorvel);

                % Channel estimation for each anchor
                X = cell(1, numAnchors);
                for idxAnchor = 1:numAnchors
                    bpskSymbol = randi([0,1], [N M]) * 2 - 1;

                    % OFDM modulation (without Communications Toolbox)
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

                % --- Estimation ---
                if strcmp(measurementType, 'TOA')
                    if strcmp(spectrumMethod, 'FFT')
                        estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                            'Measurement', 'TOA', 'SpectrumMethod', 'FFT', ...
                            'VarianceOutputPort', true, 'DelayOffsetInputPort', true);
                    else
                        estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                            'Measurement', 'TOA', 'SpectrumMethod', 'MUSIC', ...
                            'VarianceOutputPort', true, 'DelayOffsetInputPort', true, ...
                            'ForwardBackwardAveraging', true, 'SpatialSmoothing', ceil(N/2));
                    end
                    [Y, estVar] = estimator(X, freqSpacing, delayoffset);
                    tgtposEstAll(:, idxTgt) = toaposest(Y, estVar, anchorpos);
                else
                    if strcmp(spectrumMethod, 'FFT')
                        estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                            'Measurement', 'TDOA', 'SpectrumMethod', 'FFT', ...
                            'VarianceOutputPort', true, 'DelayOffsetInputPort', true);
                    else
                        estimator = phased.TOAEstimator('PropagationSpeed', cLight, ...
                            'Measurement', 'TDOA', 'SpectrumMethod', 'MUSIC', ...
                            'VarianceOutputPort', true, 'DelayOffsetInputPort', true, ...
                            'ForwardBackwardAveraging', true, 'SpatialSmoothing', ceil(N/2));
                    end
                    [Y, estVar] = estimator(X, freqSpacing, delayoffset);
                    tgtposEstAll(:, idxTgt) = tdoaposest(Y, estVar, anchorpos);
                end
            end

            % --- Compute RMSE ---
            errors = sqrt(sum((tgtposEstAll - tgtposAll).^2, 1));
            rmseVal = sqrt(mean(errors.^2));

            % --- Update RMSE label ---
            lblRMSE.Text = sprintf('RMSE: %.4f m', rmseVal);

            % --- Plot results ---
            cla(ax);
            hold(ax, 'on');

            % Anchors
            plot(ax, anchorpos(1,:), anchorpos(2,:), 'b^', ...
                'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Anchors');

            % True target positions
            plot(ax, tgtposAll(1,:), tgtposAll(2,:), 'rx', ...
                'LineWidth', 2, 'MarkerSize', 12, 'DisplayName', 'True Positions');

            % Estimated positions
            plot(ax, tgtposEstAll(1,:), tgtposEstAll(2,:), 'go', ...
                'LineWidth', 2, 'MarkerSize', 12, 'DisplayName', 'Estimated Positions');

            % Lines connecting true -> estimated
            for k = 1:numTargets
                plot(ax, [tgtposAll(1,k), tgtposEstAll(1,k)], ...
                         [tgtposAll(2,k), tgtposEstAll(2,k)], ...
                    'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            end

            legend(ax, 'Location', 'best', 'FontSize', 10);
            title(ax, sprintf('%s Localization (%s) — RMSE: %.4f m', ...
                measurementType, spectrumMethod, rmseVal));
            xlabel(ax, 'X (meters)');
            ylabel(ax, 'Y (meters)');
            grid(ax, 'on');
            hold(ax, 'off');

            lblStatus.Text = sprintf('Done. %s/%s — %d anchors, %d target(s).', ...
                measurementType, spectrumMethod, numAnchors, numTargets);

        catch ME
            lblStatus.Text = ['Error: ', ME.message];
            lblRMSE.Text = 'RMSE: Error';
        end

        % --- Re-enable Run button ---
        btnRun.Enable = 'on';
        btnRun.Text = 'Run Simulation';
    end
end