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

    % --- Show Target Curves ---
    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Curves:', ...
        'Position', [leftPad yPos 55 ctrlHeight], ...
        'FontWeight', 'bold');
    chkTargets = gobjects(1, 5);
    chkW = 45;
    for k = 1:5
        kk = k;
        chkTargets(k) = uicheckbox(controlPanel, ...
            'Text', sprintf('T%d', k), ...
            'Value', true, ...
            'Visible', 'off', ...
            'Position', [leftPad + 55 + (k-1)*chkW, yPos, chkW, ctrlHeight], ...
            'ValueChangedFcn', @(~,~) toggleTarget(kk));
    end

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
            YAll = cell(1, numTargets);
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
                YAll{idxTgt} = Y;
            end

            % --- Compute RMSE ---
            errors = sqrt(sum((tgtposEstAll - tgtposAll).^2, 1));
            rmseVal = sqrt(mean(errors.^2));

            % --- Update RMSE label ---
            lblRMSE.Text = sprintf('RMSE: %.4f m', rmseVal);

            % --- Plot results ---
            delete(allchild(ax));
            hold(ax, 'on');

            % Anchors
            plot(ax, anchorpos(1,:), anchorpos(2,:), 'b^', ...
                'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Anchors');

            % Per-target colors (shared by markers, connecting lines, and curves)
            tgtColors = lines(numTargets);

            % True target positions, estimated positions, connecting lines, labels
            for k = 1:numTargets
                tgtColor = tgtColors(k, :);
                tgtTag = sprintf('Target%d', k);
                if k == 1
                    plot(ax, tgtposAll(1,k), tgtposAll(2,k), 'x', ...
                        'Color', tgtColor, 'LineWidth', 2, 'MarkerSize', 12, ...
                        'Tag', tgtTag, 'DisplayName', 'True Positions');
                    plot(ax, tgtposEstAll(1,k), tgtposEstAll(2,k), 'o', ...
                        'Color', tgtColor, 'LineWidth', 2, 'MarkerSize', 12, ...
                        'Tag', tgtTag, 'DisplayName', 'Estimated Positions');
                else
                    plot(ax, tgtposAll(1,k), tgtposAll(2,k), 'x', ...
                        'Color', tgtColor, 'LineWidth', 2, 'MarkerSize', 12, ...
                        'Tag', tgtTag, 'HandleVisibility', 'off');
                    plot(ax, tgtposEstAll(1,k), tgtposEstAll(2,k), 'o', ...
                        'Color', tgtColor, 'LineWidth', 2, 'MarkerSize', 12, ...
                        'Tag', tgtTag, 'HandleVisibility', 'off');
                end
                plot(ax, [tgtposAll(1,k), tgtposEstAll(1,k)], ...
                         [tgtposAll(2,k), tgtposEstAll(2,k)], ...
                    '--', 'Color', tgtColor, 'LineWidth', 1, ...
                    'Tag', tgtTag, 'HandleVisibility', 'off');
                % Target number label (offset slightly above the true position)
                text(ax, tgtposAll(1,k), tgtposAll(2,k), sprintf('  T%d', k), ...
                    'Color', tgtColor, 'FontWeight', 'bold', 'FontSize', 10, ...
                    'Tag', tgtTag, 'VerticalAlignment', 'bottom');
            end

            % Trilateration circles (TOA) or hyperbola curves (TDOA)
            angles = 0:2*pi/720:2*pi;
            if strcmp(measurementType, 'TOA')
                curveName = 'Trilateration Circles';
            else
                curveName = 'Hyperbola Curves';
            end
            for idxTgt = 1:numTargets
                tgtColor = tgtColors(idxTgt, :);
                curveTag = sprintf('Curves_Target%d', idxTgt);
                if strcmp(measurementType, 'TOA')
                    rngEst = YAll{idxTgt} * cLight;
                    for anchorIdx = 1:numAnchors
                        cx = rngEst(anchorIdx) * cos(angles) + anchorpos(1, anchorIdx);
                        cy = rngEst(anchorIdx) * sin(angles) + anchorpos(2, anchorIdx);
                        if anchorIdx == 1
                            plot(ax, cx, cy, '--', 'Color', tgtColor, 'LineWidth', 1, ...
                                'Tag', curveTag, ...
                                'DisplayName', sprintf('%s (Target %d)', curveName, idxTgt));
                        else
                            plot(ax, cx, cy, '--', 'Color', tgtColor, 'LineWidth', 1, ...
                                'Tag', curveTag, ...
                                'HandleVisibility', 'off');
                        end
                    end
                else
                    rngDiffEst = YAll{idxTgt} * cLight;
                    numAnchorPair = length(rngDiffEst);
                    firstCurveForTgt = true;
                    for anchorPairIdx = 1:numAnchorPair
                        [hx, hy] = get2DHyperbolicSurface( ...
                            anchorpos(:, 1), anchorpos(:, anchorPairIdx+1), rngDiffEst(anchorPairIdx));
                        if isreal(hx) && isreal(hy)
                            if firstCurveForTgt
                                plot(ax, hx, hy, '--', 'Color', tgtColor, 'LineWidth', 1, ...
                                    'Tag', curveTag, ...
                                    'DisplayName', sprintf('%s (Target %d)', curveName, idxTgt));
                                firstCurveForTgt = false;
                            else
                                plot(ax, hx, hy, '--', 'Color', tgtColor, 'LineWidth', 1, ...
                                    'Tag', curveTag, ...
                                    'HandleVisibility', 'off');
                            end
                        end
                    end
                end
            end

            % Refresh per-target checkboxes
            for k = 1:5
                if k <= numTargets
                    chkTargets(k).Value = true;
                    chkTargets(k).Visible = 'on';
                else
                    chkTargets(k).Visible = 'off';
                end
            end

            legend(ax, 'Location', 'best', 'FontSize', 10);
            title(ax, sprintf('%s Localization (%s) — RMSE: %.4f m', ...
                measurementType, spectrumMethod, rmseVal));
            xlabel(ax, 'X (meters)');
            ylabel(ax, 'Y (meters)');
            grid(ax, 'on');

            % Set axis limits based on anchors and targets (avoids hyperbola infinities)
            allX = [anchorpos(1,:), tgtposAll(1,:), tgtposEstAll(1,:)];
            allY = [anchorpos(2,:), tgtposAll(2,:), tgtposEstAll(2,:)];
            pad = max([max(allX)-min(allX), max(allY)-min(allY)]) * 0.3;
            pad = max(pad, 5);
            xlim(ax, [min(allX) - pad, max(allX) + pad]);
            ylim(ax, [min(allY) - pad, max(allY) + pad]);

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

    function toggleTarget(tgtIdx)
        objs = findall(ax, 'Tag', sprintf('Curves_Target%d', tgtIdx));
        if isempty(objs)
            return;
        end
        if chkTargets(tgtIdx).Value
            set(objs, 'Visible', 'on');
        else
            set(objs, 'Visible', 'off');
        end
    end
end

function [x, y] = get2DHyperbolicSurface(anchorRefPos, anchorPos, rngDiffEst)
% Get 2D hyperbolic surface for a given pair of anchors
theta = linspace(-pi/2 * 0.98, pi/2 * 0.98, 300);
phi = 0;
[Theta, Phi] = meshgrid(theta, phi);

D = norm(anchorRefPos - anchorPos) / 2;
c = rngDiffEst;

xC = -c ./ cos(Theta) ./ 2;
yC = sqrt(4*D^2 - c^2) .* tan(Theta) .* cos(Phi) ./ 2;
zC = sqrt(4*D^2 - c^2) .* tan(Theta) .* sin(Phi) ./ 2;

r0 = (anchorPos + anchorRefPos) / 2;
a = [1; 0; 0];
b = (anchorPos - anchorRefPos);
b = b / norm(b);
v = cross(a, b);
s = norm(v);
c = dot(a, b);
V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
if abs(s) > 0
    R = eye(3) + V + V^2 * (1 - c) / s^2;
else
    R = eye(3);
end

x = R(1,1).*xC + R(1,2).*yC + R(1,3).*zC;
y = R(2,1).*xC + R(2,2).*yC + R(2,3).*zC;
x = x + r0(1);
y = y + r0(2);
end