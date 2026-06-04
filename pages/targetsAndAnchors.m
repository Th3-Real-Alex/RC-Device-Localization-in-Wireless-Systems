function targetsAndAnchors()
    pageContainer = uifigure('Name', 'Device Localization in Wireless Systems - Target & Anchor Effects', ...
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
    rowSpacing = 32;
    halfW = floor(ctrlWidth/2) - 5;

    % ---- Left Panel (controls) ----
    controlPanel = uipanel(pageContainer, ...
        'Title', 'Sweep Parameters', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'Position', [10 10 panelWidth 700]);

    yPos = 645;

    yPos = yPos - 5;
    uilabel(controlPanel, 'Text', 'Measurement Type', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddMeasurement = uidropdown(controlPanel, ...
        'Items', {'TOA', 'TDOA', 'Both'}, 'Value', 'Both', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Spectrum Method', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    ddSpectrum = uidropdown(controlPanel, ...
        'Items', {'FFT', 'MUSIC'}, 'Value', 'FFT', ...
        'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Anchor sweep min / max (3-10)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnAnchorMin = uispinner(controlPanel, 'Value', 3, 'Limits', [3 10], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos halfW ctrlHeight]);
    spnAnchorMax = uispinner(controlPanel, 'Value', 10, 'Limits', [3 10], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad+halfW+10 yPos halfW ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Targets while sweeping anchors', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnFixedTargets = uispinner(controlPanel, 'Value', 1, 'Limits', [1 5], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos ctrlWidth ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Target sweep min / max (1-5)', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnTargetMin = uispinner(controlPanel, 'Value', 1, 'Limits', [1 5], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos halfW ctrlHeight]);
    spnTargetMax = uispinner(controlPanel, 'Value', 5, 'Limits', [1 5], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad+halfW+10 yPos halfW ctrlHeight]);

    yPos = yPos - rowSpacing;
    uilabel(controlPanel, 'Text', 'Anchors while sweeping targets', ...
        'Position', [leftPad yPos ctrlWidth labelHeight], 'FontWeight', 'bold');
    yPos = yPos - ctrlHeight - 2;
    spnFixedAnchors = uispinner(controlPanel, 'Value', 5, 'Limits', [3 10], 'Step', 1, ...
        'RoundFractionalValues', 'on', 'Position', [leftPad yPos ctrlWidth ctrlHeight]);

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

    yPos = yPos - 70;
    lblSummary = uilabel(controlPanel, ...
        'Text', 'Summary: —', ...
        'FontSize', 11, 'FontWeight', 'bold', 'FontColor', [0.2 0.2 0.8], ...
        'WordWrap', 'on', ...
        'Position', [leftPad yPos ctrlWidth 65]);

    % ---- Right Panel (two stacked axes) ----
    axW = 1280 - panelWidth - 60;
    axH = 295;
    axAnchors = uiaxes(pageContainer, 'Position', [panelWidth + 30, 380, axW, axH]);
    axTargets = uiaxes(pageContainer, 'Position', [panelWidth + 30, 40, axW, axH]);
    title(axAnchors, 'RMSE vs Number of Anchors');
    xlabel(axAnchors, 'Number of Anchors'); ylabel(axAnchors, 'RMSE (m)');
    grid(axAnchors, 'on');
    title(axTargets, 'RMSE vs Number of Targets');
    xlabel(axTargets, 'Number of Targets'); ylabel(axTargets, 'RMSE (m)');
    grid(axTargets, 'on');

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
            measSel = ddMeasurement.Value;
            if strcmp(measSel, 'Both')
                measList = {'TOA', 'TDOA'};
            else
                measList = {measSel};
            end
            spectrumMethod = ddSpectrum.Value;
            anchorMin = spnAnchorMin.Value;
            anchorMax = max(spnAnchorMax.Value, anchorMin);
            targetMin = spnTargetMin.Value;
            targetMax = max(spnTargetMax.Value, targetMin);

            baseParams = struct();
            baseParams.spectrumMethod = spectrumMethod;
            baseParams.bw          = spnBandwidth.Value * 1e6;
            baseParams.Pt          = spnPower.Value;
            baseParams.NF          = spnNoise.Value;
            baseParams.delayoffset = 0;
            baseParams.N           = str2double(ddSubcarriers.Value);
            baseParams.M           = str2double(ddSymbols.Value);

            fixedTargets = spnFixedTargets.Value;
            fixedAnchors = spnFixedAnchors.Value;
            numTrials    = spnTrials.Value;

            anchorSweep = anchorMin:anchorMax;
            targetSweep = targetMin:targetMax;

            if strcmp(spectrumMethod, 'MUSIC')
                lblStatus.Text = 'MUSIC selected — sweep may take several minutes...';
                drawnow;
            end

            cla(axAnchors); hold(axAnchors, 'on'); grid(axAnchors, 'on');
            cla(axTargets); hold(axTargets, 'on'); grid(axTargets, 'on');

            colorMap = struct('TOA', [0.00 0.45 0.74], 'TDOA', [0.85 0.33 0.10]);

            % ----- Sweep over anchors -----
            anchorMeans = nan(numel(measList), numel(anchorSweep));
            for mi = 1:numel(measList)
                meas = measList{mi};
                minAnchorsForMeas = 3 + strcmp(meas, 'TDOA');
                rmseMat = nan(numTrials, numel(anchorSweep));
                for i = 1:numel(anchorSweep)
                    nA = anchorSweep(i);
                    if nA < minAnchorsForMeas
                        continue;
                    end
                    for t = 1:numTrials
                        lblStatus.Text = sprintf('[%s] Anchors %d/%d  trial %d/%d', ...
                            meas, i, numel(anchorSweep), t, numTrials);
                        drawnow;
                        p = baseParams;
                        p.numAnchors      = nA;
                        p.numTargets      = fixedTargets;
                        p.measurementType = meas;
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
                anchorMeans(mi, :) = meanRmse;

                col = colorMap.(meas);
                fill(axAnchors, [anchorSweep fliplr(anchorSweep)], ...
                    [minRmse fliplr(maxRmse)], col, ...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(axAnchors, anchorSweep, meanRmse, '-o', ...
                    'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, ...
                    'DisplayName', sprintf('%s mean RMSE', meas));
            end
            legend(axAnchors, 'Location', 'best');

            % ----- Sweep over targets -----
            targetMeans = nan(numel(measList), numel(targetSweep));
            for mi = 1:numel(measList)
                meas = measList{mi};
                rmseMat = nan(numTrials, numel(targetSweep));
                for i = 1:numel(targetSweep)
                    nT = targetSweep(i);
                    for t = 1:numTrials
                        lblStatus.Text = sprintf('[%s] Targets %d/%d  trial %d/%d', ...
                            meas, i, numel(targetSweep), t, numTrials);
                        drawnow;
                        p = baseParams;
                        p.numAnchors      = fixedAnchors;
                        p.numTargets      = nT;
                        p.measurementType = meas;
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
                targetMeans(mi, :) = meanRmse;

                col = colorMap.(meas);
                fill(axTargets, [targetSweep fliplr(targetSweep)], ...
                    [minRmse fliplr(maxRmse)], col, ...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(axTargets, targetSweep, meanRmse, '-o', ...
                    'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, ...
                    'DisplayName', sprintf('%s mean RMSE', meas));
            end
            legend(axTargets, 'Location', 'best');

            hold(axAnchors, 'off'); hold(axTargets, 'off');

            summaryLines = strings(numel(measList), 1);
            for mi = 1:numel(measList)
                aRow = anchorMeans(mi,:);
                tRow = targetMeans(mi,:);
                aFirst = find(isfinite(aRow), 1, 'first');
                aLast  = find(isfinite(aRow), 1, 'last');
                tFirst = find(isfinite(tRow), 1, 'first');
                tLast  = find(isfinite(tRow), 1, 'last');
                if isempty(aFirst) || isempty(tFirst)
                    summaryLines(mi) = sprintf('%s: no valid results', measList{mi});
                else
                    summaryLines(mi) = sprintf('%s: anchors %.2f→%.2f m, targets %.2f→%.2f m', ...
                        measList{mi}, aRow(aFirst), aRow(aLast), tRow(tFirst), tRow(tLast));
                end
            end
            lblSummary.Text = char(strjoin(["Summary:"; summaryLines], newline));
            lblStatus.Text = sprintf('Done. %s, %s, %d trials/point.', ...
                strjoin(measList, '+'), spectrumMethod, numTrials);

        catch ME
            lblStatus.Text = ['Error: ', ME.message];
            lblSummary.Text = 'Summary: error';
        end

        btnRun.Enable = 'on';
        btnRun.Text = 'Run Sweep';
    end
end
