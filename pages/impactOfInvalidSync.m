function impactOfInvalidSync()
    pageContainer = uifigure('Name', 'Device Localization in Wireless Systems - Impact of invalid synchronization between anchors', ...
        'NumberTitle', 'off', ...
        'Resize', 'off', ...
        'ToolBar', 'none', ...
        'Position', [200 200 1280 720]);

    % ---- Constants ----
    panelWidth = 310;
    leftPad = 15;
    ctrlWidth = 280;
    ctrlHeight = 22;
    labelHeight = 16;
    rowSpacing = 34;
    halfW = floor(ctrlWidth/2) - 5;
    third = floor(ctrlWidth/3) - 6;

    controlPanel = uipanel(pageContainer, ...
        'Title', 'Sweep Parameters', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'Position', [10 10 panelWidth 700]);

    yPos = 645;

    yPos = yPos - 5;
    uilabel(controlPanel, 'Text', 'Spectrum Method', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSpectrum = uidropdown(controlPanel, ...
        'Items', {'FFT', 'MUSIC'}, 'Value', 'FFT', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Anchors  /  Targets', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnAnchors = uispinner(controlPanel, 'Value', 5, 'Limits', [3 10], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos halfW ctrlHeight]);
    spnTargets = uispinner(controlPanel, 'Value', 1, 'Limits', [1 5], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad+halfW+10 yPos halfW ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Delay offset min / max / step (ns)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnDelayMin = uispinner(controlPanel, 'Value', 0, 'Limits', [0 1000], 'Step', 10, ...
        'Position', [leftPad yPos third ctrlHeight]);
    spnDelayMax = uispinner(controlPanel, 'Value', 500, 'Limits', [0 1000], 'Step', 10, ...
        'Position', [leftPad+third+5 yPos third ctrlHeight]);
    spnDelayStep = uispinner(controlPanel, 'Value', 50, 'Limits', [1 500], 'Step', 10, ...
        'Position', [leftPad+2*(third+5) yPos third ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Bandwidth (MHz)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnBandwidth = uispinner(controlPanel, 'Value', 100, 'Limits', [50 400], 'Step', 50, ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Tx Power (W)  /  Noise Figure (dB)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnPower = uispinner(controlPanel, 'Value', 0.1, 'Limits', [0.01 1], 'Step', 0.05, ...
        'Position', [leftPad yPos halfW ctrlHeight]);
    spnNoise = uispinner(controlPanel, 'Value', 2.9, 'Limits', [1 10], 'Step', 0.5, ...
        'Position', [leftPad+halfW+10 yPos halfW ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Subcarriers (N)  /  OFDM Symbols (M)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSubcarriers = uidropdown(controlPanel, 'Items', {'256','512','1024','2048'}, ...
        'Value', '1024', 'Position', [leftPad yPos halfW ctrlHeight]);
    ddSymbols = uidropdown(controlPanel, 'Items', {'4','8','16'}, ...
        'Value', '8', 'Position', [leftPad+halfW+10 yPos halfW ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Monte Carlo Trials per Point', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnTrials = uispinner(controlPanel, 'Value', 10, 'Limits', [1 30], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing - 5;
    btnRun = uibutton(controlPanel, ...
        'Text', 'Run Sweep', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.93 0.53 0.24], 'FontColor', 'white', ...
        'Position', [leftPad yPos ctrlWidth 35], ...
        'ButtonPushedFcn', @(~,~) runSweep());

    yPos = yPos - 95;
    uilabel(controlPanel, ...
        'Text', ['Expected: TOA error grows ≈ 0.3 m / ns (range bias = c·Δt). ', ...
                 'TDOA cancels a common offset and stays roughly flat.'], ...
        'FontSize', 11, 'FontColor', [0.3 0.3 0.3], 'WordWrap', 'on', ...
        'Position', [leftPad yPos ctrlWidth 60]);

    yPos = yPos - 60;
    lblSummary = uilabel(controlPanel, ...
        'Text', 'Summary: —', ...
        'FontSize', 11, 'FontWeight', 'bold', 'FontColor', [0.2 0.2 0.8], ...
        'WordWrap', 'on', ...
        'Position', [leftPad yPos ctrlWidth 55]);

    % ---- Right Panel ----
    ax = uiaxes(pageContainer, ...
        'Position', [panelWidth + 30, 40, 1280 - panelWidth - 60, 640]);
    title(ax, 'RMSE vs Synchronization Offset');
    xlabel(ax, 'Delay offset (ns)');
    ylabel(ax, 'RMSE (m)');
    grid(ax, 'on');

    lblStatus = uilabel(pageContainer, ...
        'Text', 'Ready. Configure parameters and click Run Sweep.', ...
        'FontSize', 12, 'FontColor', [0.3 0.3 0.3], ...
        'Position', [panelWidth + 30, 690, 900, 25]);

    % ======================================================================
    function runSweep()
        btnRun.Enable = 'off';
        btnRun.Text = 'Running...';
        lblSummary.Text = 'Summary: —';
        drawnow;

        try
            spectrumMethod = ddSpectrum.Value;
            numAnchors  = spnAnchors.Value;
            numTargets  = spnTargets.Value;
            delayMin    = spnDelayMin.Value;
            delayMax    = max(spnDelayMax.Value, delayMin);
            delayStep   = max(spnDelayStep.Value, 1);
            numTrials   = spnTrials.Value;

            baseParams = struct();
            baseParams.spectrumMethod = spectrumMethod;
            baseParams.bw          = spnBandwidth.Value * 1e6;
            baseParams.Pt          = spnPower.Value;
            baseParams.NF          = spnNoise.Value;
            baseParams.N           = str2double(ddSubcarriers.Value);
            baseParams.M           = str2double(ddSymbols.Value);
            baseParams.numAnchors  = numAnchors;
            baseParams.numTargets  = numTargets;

            offsetsNs = delayMin:delayStep:delayMax;
            if isempty(offsetsNs)
                offsetsNs = delayMin;
            end

            if strcmp(spectrumMethod, 'MUSIC')
                lblStatus.Text = 'MUSIC selected — sweep may take several minutes...';
                drawnow;
            end

            cla(ax); hold(ax, 'on'); grid(ax, 'on');

            measList = {'TOA', 'TDOA'};
            colorMap = struct('TOA', [0.00 0.45 0.74], 'TDOA', [0.85 0.33 0.10]);
            meanCurves = nan(numel(measList), numel(offsetsNs));

            for mi = 1:numel(measList)
                meas = measList{mi};
                rmseMat = nan(numTrials, numel(offsetsNs));
                for i = 1:numel(offsetsNs)
                    offNs = offsetsNs(i);
                    for t = 1:numTrials
                        lblStatus.Text = sprintf('[%s] offset %g ns (%d/%d)  trial %d/%d', ...
                            meas, offNs, i, numel(offsetsNs), t, numTrials);
                        drawnow;
                        p = baseParams;
                        p.measurementType = meas;
                        p.delayoffset     = offNs * 1e-9;
                        p.seed            = t;
                        try
                            rmseMat(t, i) = runLocalizationScenario(p);
                        catch
                            rmseMat(t, i) = NaN;
                        end
                    end
                end
                meanRmse = mean(rmseMat, 1, 'omitnan');
                minRmse  = min(rmseMat, [], 1, 'omitnan');
                maxRmse  = max(rmseMat, [], 1, 'omitnan');
                meanCurves(mi, :) = meanRmse;

                col = colorMap.(meas);
                if numel(offsetsNs) > 1
                    fill(ax, [offsetsNs fliplr(offsetsNs)], ...
                        [minRmse fliplr(maxRmse)], col, ...
                        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
                plot(ax, offsetsNs, meanRmse, '-o', ...
                    'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, ...
                    'DisplayName', sprintf('%s mean RMSE', meas));
            end
            legend(ax, 'Location', 'best');
            hold(ax, 'off');

            toaSlope = NaN;
            if numel(offsetsNs) > 1
                dx = offsetsNs(end) - offsetsNs(1);
                if dx > 0
                    toaSlope = (meanCurves(1, end) - meanCurves(1, 1)) / dx;
                end
            end
            lblSummary.Text = sprintf( ...
                ['Summary:\nTOA: %.2f m → %.2f m (slope ≈ %.3f m/ns)\n', ...
                 'TDOA: %.2f m → %.2f m'], ...
                meanCurves(1,1), meanCurves(1,end), toaSlope, ...
                meanCurves(2,1), meanCurves(2,end));
            lblStatus.Text = sprintf('Done. %s, %d points × %d trials.', ...
                spectrumMethod, numel(offsetsNs), numTrials);

        catch ME
            lblStatus.Text = ['Error: ', ME.message];
            lblSummary.Text = 'Summary: error';
        end

        btnRun.Enable = 'on';
        btnRun.Text = 'Run Sweep';
    end
end
