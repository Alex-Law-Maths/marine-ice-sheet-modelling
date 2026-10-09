%% ICE-SHEET MODEL DRIVER
% This driver contains the numerical routines used for the shear-dominated
% and extensional ice-sheet models. The four switches below control which
% numerical studies are run. Set a switch to false to skip that study.
%
% All solver, grounding-line, velocity-profile and plotting functions are
% included as local functions in this file; no additional MATLAB functions
% are required.

clear;
close all;
clc;

%% OUTPUT FOLDER
% Save all generated plots and parameter tables in a results folder next
% to this driver, regardless of MATLAB's current working directory.
resultsDirectory = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~isfolder(resultsDirectory)
    mkdir(resultsDirectory);
end

%% DRIVER OPTIONS
runParameterVariationProfiles   = true;
runConstantIceSheetDepthProfiles = true;
runConstantIceShelfDepthProfiles = true;
runGroundingLineParameterSweeps = true;

if runParameterVariationProfiles
    runParameterVariationProfileStudy();
end

if runConstantIceSheetDepthProfiles
    runConstantIceSheetDepthStudy();
end

if runConstantIceShelfDepthProfiles
    runConstantIceShelfDepthStudy();
end

if runGroundingLineParameterSweeps
    runGroundingLineSweepStudy();
end

%% DRIVER ROUTINES
function runParameterVariationProfileStudy()
    %% Parameter-variation profiles
    % Compare shear-dominated and extensional ice-sheet profiles for n = 3.
    % Bed slope beta, accumulation rate and slip length Lambda are varied.
    % Geometry and velocity profiles are obtained from the functions below.
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    %% COMMON PARAMETERS
    epsilon = 0.1;
    shelfExtent = 3;
    % The solver functions have q_0 = 1 built into their equations.
    qUpstream = 1;
    % true  -> plot x/x_G
    % false -> plot dimensionless x directly
    scaleXByGroundingLine = true;

    %% FIGURE EXPORT SETTINGS
    %  The exported figure dimensions may be adjusted here.
    exportSettings.figureWidthCm = 16;
    exportSettings.figureHeightCm = 8.5;

    %% MODEL ROWS
    %  Bottom: Shear-dominated ice-sheet model, Lambda = 0
    %  Middle: Shear-dominated ice-sheet model, Lambda = 0.1
    %  Top:    Extensional ice-sheet model, Lambda = 2
    % extensionalReferenceLambda
    shearDominatedLambdaRows = [0, 0.1];
    extensionalReferenceLambda = 2;
    accumulationShearDominatedLambdaRows = [0, 0.3];
    accumulationExtensionalLambda = 2;

    %% PLOT COLOURS
    colours.ice      = [0.78, 0.91, 1.00];
    colours.bedrock  = [0.48, 0.30, 0.18];
    colours.ocean    = [0.25, 0.62, 0.85];
    colours.surface  = [0.05, 0.25, 0.65];
    colours.bed      = [0.25, 0.25, 0.25];
    colours.sea      = [0.10, 0.35, 0.55];
    colours.velocity = [0.05, 0.05, 0.05];

    %% VELOCITY-ARROW SETTINGS
    velocityArrowSettings.nVerticalLevels = 3;
    velocityArrowSettings.sheetArrowXFraction = 0.45;
    velocityArrowSettings.shelfArrowXFraction = 0.45;
    velocityArrowSettings.maxLengthFraction = 0.2;
    velocityArrowSettings.colour = colours.velocity;
    velocityArrowSettings.lineWidth = 0.6;
    velocityArrowSettings.headSize = 5;

    %% BED-SLOPE VARIATION, n = 3
    betaProfileValues = [0.5, 1, 3];
    betaProfileAccumulation = 0.2;
    betaProfileN = 3;
    betaProfileResults = makeResultGrid( betaProfileValues, ...
        betaProfileAccumulation*ones(size(betaProfileValues)), ...
        betaProfileN*ones(size(betaProfileValues)), epsilon, shelfExtent, ...
        shearDominatedLambdaRows, extensionalReferenceLambda, ...
        scaleXByGroundingLine);
    plotProfileGrid( betaProfileResults, colours, velocityArrowSettings, ...
        '$\beta$', betaProfileValues, ...
        '$\Lambda$', [shearDominatedLambdaRows, extensionalReferenceLambda], ...
        scaleXByGroundingLine, exportSettings, 'bed_slope_profile_comparison.png');

    %% ACCUMULATION-RATE VARIATION, n = 3
    accumulationProfileBeta = 5;
    accumulationProfileValues = [0, 1, 5];
    accumulationProfileN = 3;
    accumulationProfileResults = makeResultGrid( ...
        accumulationProfileBeta*ones(size(accumulationProfileValues)), ...
        accumulationProfileValues, ...
        accumulationProfileN*ones(size(accumulationProfileValues)), ...
        epsilon, shelfExtent, ...
        accumulationShearDominatedLambdaRows, accumulationExtensionalLambda, ...
        scaleXByGroundingLine);
    plotProfileGrid( accumulationProfileResults, colours, velocityArrowSettings, ...
        '$\mathcal{A}$', accumulationProfileValues, '$\Lambda$', ...
        [accumulationShearDominatedLambdaRows, accumulationExtensionalLambda], ...
        scaleXByGroundingLine, exportSettings, 'accumulation_profile_comparison.png');

    %% PROFILE-PLOTTING HELPERS

    function results = makeResultGrid( ...
        betaValues, accumulationValues, nValues, epsilon, shelfExtent, ...
        shearDominatedLambdaRows, extensionalReferenceLambda, ...
        scaleXByGroundingLine)
        nColumns = numel(betaValues);
        nRows = numel(shearDominatedLambdaRows) + 1;
        results = cell(nRows, nColumns);
        for column = 1:nColumns
            beta = betaValues(column);
            accumulationRate = accumulationValues(column);
            n = nValues(column);
            for row = 1:nRows
                if row <= numel(shearDominatedLambdaRows)
                    Lambda = shearDominatedLambdaRows(row);
                    results{row,column} = callShearDominatedSolver( ...
                        Lambda, beta, accumulationRate, n, epsilon, ...
                        shelfExtent, scaleXByGroundingLine);
                else
                    results{row,column} = callExtensionalSolver( ...
                        extensionalReferenceLambda, beta, accumulationRate, n, ...
                        epsilon, shelfExtent, scaleXByGroundingLine);
                end
            end
        end
    end

    function result = callShearDominatedSolver( ...
        Lambda, beta, accumulationRate, n, epsilon, ...
        shelfExtent, scaleXByGroundingLine)
        if n == 1
            [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, ...
                hxSheet,xG,HG,hG,bG] = Sheardominatedicesheetnewtonian( ...
                    Lambda,beta,accumulationRate,epsilon,shelfExtent);
        else
            [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, ...
                hxSheet,xG,HG,hG,bG] = Sheardominatedicesheet( ...
                    Lambda,beta,accumulationRate,n,epsilon,shelfExtent);
        end
        % Store profile outputs as column vectors for interpolation and
        % plotting.
        xSheet  = xSheet(:);
        hSheet  = hSheet(:);
        bSheet  = bSheet(:);
        hxSheet = hxSheet(:);
        xShelf  = xShelf(:);
        hShelf  = hShelf(:);
        bShelf  = bShelf(:);
        % The shear-dominated ice-sheet model solver integrates the grounded
        % sheet from x_G back to 0.
        % Reverse it so interpolation and plotting use increasing x.
        if numel(xSheet) > 1 && xSheet(1) > xSheet(end)
            xSheet  = flipud(xSheet);
            hSheet  = flipud(hSheet);
            bSheet  = flipud(bSheet);
            hxSheet = flipud(hxSheet);
        end
        HSheet = hSheet - bSheet;
        HShelf = hShelf - bShelf;
        result = commonResultStruct( ...
            'shear-dominated', 'Shear-dominated ice-sheet model', ...
            beta, epsilon, accumulationRate, Lambda, n, xG, HG, hG, bG, ...
            xSheet, HSheet, hSheet, bSheet, ...
            xShelf, HShelf, hShelf, bShelf, scaleXByGroundingLine);
        result.hxSheet = hxSheet;
    end

    function result = callExtensionalSolver( ...
        Lambda, beta, accumulationRate, n, epsilon, ...
        shelfExtent, scaleXByGroundingLine)
        [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, xG,HG,hG,bG] = ...
            Extensionalicesheet( ...
                Lambda,beta,accumulationRate,n,epsilon,shelfExtent);
        % Store profile outputs as column vectors for interpolation and
        % plotting.
        xSheet = xSheet(:);
        hSheet = hSheet(:);
        bSheet = bSheet(:);
        xShelf = xShelf(:);
        hShelf = hShelf(:);
        bShelf = bShelf(:);
        % The extensional grounded sheet is integrated from x_G to 0;
        % reverse the output so that x increases.
        if numel(xSheet) > 1 && xSheet(1) > xSheet(end)
            xSheet = flipud(xSheet);
            hSheet = flipud(hSheet);
            bSheet = flipud(bSheet);
        end
        HSheet = hSheet - bSheet;
        HShelf = hShelf - bShelf;
        result = commonResultStruct( ...
            'extensional', 'Extensional ice-sheet model', ...
            beta, epsilon, accumulationRate, Lambda, n, xG, HG, hG, bG, ...
            xSheet, HSheet, hSheet, bSheet, ...
            xShelf, HShelf, hShelf, bShelf, scaleXByGroundingLine);
    end

    function result = commonResultStruct( modelKey, modelLabel, ...
        beta, epsilon, accumulationRate, Lambda, n, xG, HG, hG, bG, ...
        xSheet, HSheet, hSheet, bSheet, xShelf, HShelf, hShelf, bShelf, ...
        scaleXByGroundingLine)
        if scaleXByGroundingLine
            xScale = xG;
        else
            xScale = 1;
        end
        result.modelKey = modelKey;
        result.modelLabel = modelLabel;
        result.beta = beta;
        result.epsilon = epsilon;
        result.accumulationRate = accumulationRate;
        result.Lambda = Lambda;
        result.n = n;
        result.xG = xG;
        result.HG = HG;
        result.hG = hG;
        result.bG = bG;
        result.xN = xShelf(end);
        result.xSheet = xSheet;
        result.HSheet = HSheet;
        result.hSheet = hSheet;
        result.bSheet = bSheet;
        result.xShelf = xShelf;
        result.HShelf = HShelf;
        result.hShelf = hShelf;
        result.bShelf = bShelf;
        result.xSheetPlot = xSheet/xScale;
        result.xShelfPlot = xShelf/xScale;
        result.xGPlot = xG/xScale;
        result.xNPlot = result.xN/xScale;
        result.xMaxPlot = result.xN/xScale;
    end

    function plotProfileGrid( results, colours, arrowSettings, ...
        horizontalParameterName, horizontalParameterValues, ...
        verticalParameterName, verticalParameterValues, ...
        scaleXByGroundingLine, exportSettings, outputFilename)
        [nRows, nColumns] = size(results);
        yMinimum = inf;
        yMaximum = -inf;
        for k = 1:numel(results)
            result = results{k};
            yMinimum = min([ yMinimum; result.bSheet(:); result.bShelf(:)]);
            yMaximum = max([ yMaximum; result.hSheet(:); result.hShelf(:)]);
        end
        ySpan = yMaximum - yMinimum;
        if ySpan <= 0
            ySpan = 1;
        end
        yBottom = yMinimum - 0.14*ySpan;
        yTop    = yMaximum + 0.06*ySpan;
        maximumVelocity = max(cellfun( @profileMaximumVelocity, results(:)));
        figureWidthCm = exportSettings.figureWidthCm;
        figureHeightCm = exportSettings.figureHeightCm;
        figureHandle = figure( 'Color', 'w', 'Units', 'centimeters', ...
            'Position', [2, 2, figureWidthCm, figureHeightCm]);
        layout = tiledlayout( nRows, nColumns, 'TileSpacing', 'compact', ...
            'Padding', 'compact');
        layout.OuterPosition = [0.115, 0.155, 0.845, 0.79];
        legendHandles = [];
        for plotRow = 1:nRows
            % Show increasing Lambda from bottom to top.
            row = nRows - plotRow + 1;
            for column = 1:nColumns
                result = results{row,column};
                ax = nexttile(layout);
                handles = plotIceProfile( ...
                    ax, result, yBottom, yTop, colours, ...
                    maximumVelocity, arrowSettings);
                % Keep subplot-axis text only on the bottom-left panel.
                % This reduces clutter while retaining one reference set of
                % axes.
                if plotRow == nRows && column == 1
                    if scaleXByGroundingLine
                        xlabel(ax, '$x/x_G$', 'Interpreter', 'latex', ...
                            'FontSize', 9);
                    else
                        xlabel(ax, '$x$', 'Interpreter', 'latex', ...
                            'FontSize', 9);
                    end
                    ylabel(ax, '$z$', 'Interpreter', 'latex', 'FontSize', 9);
                    % Pull the subplot labels slightly closer to the axes.
                    ax.XLabel.Units = 'normalized';
                    ax.XLabel.Position = [0.5, -0.105, 0];
                    ax.YLabel.Units = 'normalized';
                    ax.YLabel.Position = [-0.085, 0.5, 0];
                    ax.XTickLabelMode = 'auto';
                    ax.YTickLabelMode = 'auto';
                    ax.YTick = chooseIntegerTicks(yBottom, yTop, 5);
                else
                    ax.XTickLabel = [];
                    ax.YTickLabel = [];
                end
                if plotRow == 1 && column == 1
                    legendHandles = [ handles.ice, handles.bedrock, ...
                        handles.ocean, handles.velocity];
                end
            end
        end
        sharedLegend = legend( legendHandles, ...
            {'Ice', 'Bedrock', 'Ocean', 'Velocity'}, ...
            'Interpreter', 'none', 'FontSize', 9, ...
            'Orientation', 'horizontal', 'Box', 'off');
        sharedLegend.NumColumns = 4;
        sharedLegend.Layout.Tile = 'south';
        addParameterDirectionAxes( figureHandle, ...
            horizontalParameterName, horizontalParameterValues, ...
            verticalParameterName, verticalParameterValues);
        exportgraphics( figureHandle, resultsPath(outputFilename), 'Resolution', 300);
    end

    function parameterAxes = addParameterDirectionAxes( ...
        figureHandle, horizontalName, horizontalValues, ...
        verticalName, verticalValues)
        horizontalValues = horizontalValues(:).';
        verticalValues   = verticalValues(:).';
        nHorizontal = numel(horizontalValues);
        nVertical   = numel(verticalValues);
        parameterAxes = axes( 'Parent', figureHandle, ...
            'Position', [0.105, 0.165, 0.855, 0.775], 'Color', 'none', ...
            'Box', 'off', 'XAxisLocation', 'bottom', ...
            'YAxisLocation', 'left', 'XLim', [0.5, nHorizontal + 0.5], ...
            'YLim', [0.5, nVertical + 0.5], 'XTick', 1:nHorizontal, ...
            'YTick', 1:nVertical, ...
            'XTickLabel', compose('%.3g', horizontalValues), ...
            'YTickLabel', compose('%.3g', verticalValues), ...
            'TickDir', 'out', 'TickLength', [0.012, 0.012], ...
            'LineWidth', 1.10, 'FontSize', 9, 'Layer', 'top', ...
            'HitTest', 'off', 'PickableParts', 'none');
        parameterAxes.TickLabelInterpreter = 'none';
        xlabel( parameterAxes, horizontalName, 'Interpreter', 'latex', ...
            'FontWeight', 'bold', 'FontSize', 9);
        ylabel( parameterAxes, verticalName, 'Interpreter', 'latex', ...
            'FontWeight', 'bold', 'FontSize', 9);
        % Keep the large parameter-axis labels close to the subplot grid.
        parameterAxes.XLabel.Units = 'normalized';
        parameterAxes.XLabel.Position = [0.5, -0.1, 0];
        parameterAxes.YLabel.Units = 'normalized';
        parameterAxes.YLabel.Position = [-0.065, 0.5, 0];
        parameterAxes.XGrid = 'off';
        parameterAxes.YGrid = 'off';
        parameterAxesPosition = parameterAxes.Position;
        annotation( figureHandle, 'arrow', [ ...
                parameterAxesPosition(1) + parameterAxesPosition(3), ...
                parameterAxesPosition(1) + parameterAxesPosition(3) + 0.02], ...
            [ parameterAxesPosition(2), parameterAxesPosition(2)], ...
            'LineWidth',1.10, 'HeadLength',6, 'HeadWidth',6);
        annotation( figureHandle, 'arrow', [ parameterAxesPosition(1), ...
                parameterAxesPosition(1)], [ ...
                parameterAxesPosition(2) + parameterAxesPosition(4), ...
                parameterAxesPosition(2) + parameterAxesPosition(4) + 0.02], ...
            'LineWidth',1.10, 'HeadLength',6, 'HeadWidth',6);
    end

    function yTicks = chooseIntegerTicks(yMinimum, yMaximum, maxTickCount)
        if nargin < 3
            maxTickCount = 5;
        end
        startTick = ceil(yMinimum);
        endTick = floor(yMaximum);
        if startTick > endTick
            yTicks = round(linspace(yMinimum, yMaximum, min(3,maxTickCount)));
            yTicks = unique(yTicks);
            return;
        end
        integerTicks = startTick:endTick;
        if numel(integerTicks) <= maxTickCount
            yTicks = integerTicks;
            return;
        end
        step = ceil((numel(integerTicks) - 1)/(maxTickCount - 1));
        yTicks = integerTicks(1:step:end);
        if yTicks(end) ~= integerTicks(end)
            yTicks(end + 1) = integerTicks(end);
        end
    end

    function handles = plotIceProfile( ax, result, yBottom, yTop, colours, ...
        maximumVelocity, arrowSettings)
        hold(ax, 'on');
        % Treat all plotted profiles as column vectors so patch/vertcat
        % receives arrays with consistent dimensions.
        xSheet = result.xSheetPlot(:);
        xShelf = result.xShelfPlot(:);
        hSheet = result.hSheet(:);
        bSheet = result.bSheet(:);
        hShelf = result.hShelf(:);
        bShelf = result.bShelf(:);
        %% Continue bed beneath the floating shelf
        if result.xGPlot ~= 0
            xScale = result.xG/result.xGPlot;
        else
            xScale = 1;
        end
        xShelfPhysical = xShelf.*xScale;
        bOcean = -result.beta.*xShelfPhysical;
        handles.bedrock = patch( ax, [xSheet; flipud(xSheet)], ...
            [bSheet; yBottom.*ones(size(xSheet))], colours.bedrock, ...
            'EdgeColor', 'none');
        patch( ax, [xShelf; flipud(xShelf)], ...
            [bOcean; yBottom.*ones(size(xShelf))], colours.bedrock, ...
            'EdgeColor', 'none');
        handles.ocean = patch( ax, [xShelf; flipud(xShelf)], ...
            [bShelf; flipud(bOcean)], colours.ocean, 'EdgeColor', 'none');
        handles.ice = patch( ax, [xSheet; flipud(xSheet)], ...
            [hSheet; flipud(bSheet)], colours.ice, 'EdgeColor', 'none');
        patch( ax, [xShelf; flipud(xShelf)], [hShelf; flipud(bShelf)], ...
            colours.ice, 'EdgeColor', 'none');
        handles.velocity = addVelocityProfiles( ...
            ax, result, maximumVelocity, arrowSettings);
        plot(ax, xSheet, hSheet, 'Color', colours.surface, 'LineWidth', 1.5);
        plot(ax, xShelf, hShelf, 'Color', colours.surface, 'LineWidth', 1.5);
        plot(ax, xSheet, bSheet, 'Color', colours.bed, 'LineWidth', 1.5);
        plot(ax, xShelf, bShelf, 'Color', colours.surface, 'LineWidth', 1.5);
        plot(ax, xShelf, bOcean, 'Color', colours.bed, 'LineWidth', 1.5);
        xline(ax, result.xGPlot, '--k', 'LineWidth', 0.9);
        yline(ax, 0, '--', 'Color', colours.sea, 'LineWidth', 0.9);
        xlim(ax, [0, result.xMaxPlot]);
        ylim(ax, [yBottom, yTop]);
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 9);
        box(ax, 'on');
        hold(ax, 'off');
    end

    function maximumVelocity = profileMaximumVelocity(result)
        if strcmp(result.modelKey,'shear-dominated')
            etaSurface = ones(size(result.HSheet));
            uSheet = velocityprofileSheardominatedsheet( etaSurface, ...
                result.HSheet, result.hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            etaSurface = ones(size(result.HSheet));
            uSheet = velocityprofileExtensionalsheet( etaSurface, ...
                result.xSheet, result.HSheet, result.accumulationRate);
        else
            error('Unknown model key: %s.',result.modelKey);
        end
        etaShelf = ones(size(result.HShelf));
        uShelf = velocityprofileshelf( etaShelf, result.xShelf, ...
            result.HShelf, result.accumulationRate);
        maximumVelocity = max([ uSheet(:); uShelf(:)]);
    end

    function arrowHandle = addVelocityProfiles( ...
        ax, result, maximumVelocity, settings)
        if maximumVelocity <= 0 || ~isfinite(maximumVelocity)
            error('maximumVelocity must be positive and finite.');
        end
        etaValues = linspace( 0.08, 0.92, settings.nVerticalLevels).';
        if result.xGPlot ~= 0
            xScale = result.xG/result.xGPlot;
        else
            xScale = 1;
        end
        xSheetPhysical = settings.sheetArrowXFraction * result.xG;
        xSheetPosition = xSheetPhysical/xScale;
        HSheet = interp1( result.xSheet, result.HSheet, ...
            xSheetPhysical, 'pchip');
        bSheet = interp1( result.xSheet, result.bSheet, ...
            xSheetPhysical, 'pchip');
        zSheet = bSheet + etaValues.*HSheet;
        if strcmp(result.modelKey,'shear-dominated')
            hxSheet = interp1( result.xSheet, result.hxSheet, ...
                xSheetPhysical, 'pchip');
            uSheet = velocityprofileSheardominatedsheet( etaValues, ...
                HSheet, hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            uSheet = velocityprofileExtensionalsheet( etaValues, ...
                xSheetPhysical, HSheet, result.accumulationRate);
        else
            error('Unknown model key: %s.',result.modelKey);
        end
        xShelfPhysical = result.xG + settings.shelfArrowXFraction ...
              * (result.xN - result.xG);
        xShelfPosition = xShelfPhysical/xScale;
        HShelf = interp1( result.xShelf, result.HShelf, ...
            xShelfPhysical, 'pchip');
        bShelf = interp1( result.xShelf, result.bShelf, ...
            xShelfPhysical, 'pchip');
        zShelf = bShelf + etaValues.*HShelf;
        uShelf = velocityprofileshelf( etaValues, xShelfPhysical, HShelf, ...
            result.accumulationRate);
        arrowX = [ xSheetPosition.*ones(size(etaValues)); ...
            xShelfPosition.*ones(size(etaValues))];
        arrowZ = [zSheet; zShelf];
        velocityValues = [uSheet(:); uShelf(:)];
        maximumArrowLength = settings.maxLengthFraction * result.xMaxPlot;
        plottedArrowLengths = maximumArrowLength .* velocityValues ...
            ./ maximumVelocity;
        arrowHandle = quiver( ax, arrowX, arrowZ, ...
            plottedArrowLengths, zeros(size(plottedArrowLengths)), 0, ...
            'Color', settings.colour, 'LineWidth', settings.lineWidth, ...
            'MaxHeadSize', settings.headSize, 'AutoScale', 'off');
    end
end

function runConstantIceSheetDepthStudy()
    %% Constant ice-sheet-depth profiles
    % Compare constant-depth grounded ice-sheet solutions for n = 1 and n = 3.
    % Profiles are generated from the solver and velocity functions below.
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    %% MODEL PARAMETERS
    epsilon = 0.1;
    LambdaShearDominated = 0.1;
    LambdaExtensional = 2;
    qUpstream = 1;
    beta = 1;
    accumulationRate = 0;
    nValues = [1, 3];
    shelfExtent = 3;
    scaleXByGroundingLine = true;

    %% PLOT COLOURS
    colours.ice      = [0.78, 0.91, 1.00];
    colours.bedrock  = [0.48, 0.30, 0.18];
    colours.ocean    = [0.25, 0.62, 0.85];
    colours.surface  = [0.05, 0.25, 0.65];
    colours.bed      = [0.25, 0.25, 0.25];
    colours.sea      = [0.10, 0.35, 0.55];
    colours.velocity = [0.05, 0.05, 0.05];

    %% VELOCITY-ARROW SETTINGS
    velocityArrowSettings.nVerticalLevels = 3;
    velocityArrowSettings.sheetArrowXFraction = 0.45;
    velocityArrowSettings.shelfArrowXFraction = 0.45;
    velocityArrowSettings.maxLengthFraction = 0.11;
    velocityArrowSettings.colour = colours.velocity;
    velocityArrowSettings.lineWidth = 0.85;
    velocityArrowSettings.headSize = 0.75;
    exportSettings.figureWidthCm = 16;
    exportSettings.figureHeightCm = 8.0;

    %% BUILD THE FOUR PROFILES
    %  Rows:
    %    1 -> Extensional ice-sheet model, Lambda = 2
    %    2 -> Shear-dominated ice-sheet model, Lambda = 0.1
    %  Cols:
    %    1 -> n = 1
    %    2 -> n = 3
    nRows = 2;
    nColumns = numel(nValues);
    results = cell(nRows,nColumns);
    for column = 1:nColumns
        n = nValues(column);
        HGExtensionalConstant = ( ((n+1)/(n+2))^n / LambdaExtensional ...
            )^(1/(n^2 + 3*n + 1));
        betaExtensional = ( 1 ...
                / (LambdaExtensional*HGExtensionalConstant^(n+1)) )^(1/n);
        results{1,column} = callExtensionalSolver( betaExtensional, ...
            epsilon, accumulationRate, LambdaExtensional, n, shelfExtent, ...
            scaleXByGroundingLine);
        results{2,column} = callShearDominatedSolver( beta, epsilon, ...
            accumulationRate, LambdaShearDominated, n, shelfExtent, ...
            scaleXByGroundingLine);
    end

    %% EXPORT PARAMETERS USED FOR CONSTANT SHEET DEPTH
    %  Two files are written:
    %      constant_sheet_depth_parameters.csv
    %      constant_sheet_depth_parameters.tex
    %  The .tex file contains only the tabular environment, so it can be
    %  inserted directly into a LaTeX table using \input{...}.
    exportConstantSheetDepthParameters( results, ...
        'constant_sheet_depth_parameters.csv', ...
        'constant_sheet_depth_parameters.tex');

    %% PLOT
    plotConstantDepthProfilesGrid( results, colours, ...
        velocityArrowSettings, epsilon, nValues, ...
        [LambdaShearDominated, LambdaExtensional], scaleXByGroundingLine, ...
        exportSettings, 'constant_ice_sheet_depth_profiles.png');

    %% CONSTANT-DEPTH PROFILE HELPERS

    function result = callShearDominatedSolver( ...
        beta,epsilon,Accumulation,Lambda,n, shelfExtent,scaleXByGroundingLine)
        if n == 1
            [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, ...
                hxSheet,xG,HG,hG,bG] = Sheardominatedicesheetnewtonian( ...
                    Lambda, beta, Accumulation, epsilon, shelfExtent);
        else
            [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, ...
                hxSheet,xG,HG,hG,bG] = Sheardominatedicesheet( Lambda, ...
                    beta, Accumulation, n, epsilon, shelfExtent);
        end
        xSheet  = xSheet(:);
        hSheet  = hSheet(:);
        bSheet  = bSheet(:);
        hxSheet = hxSheet(:);
        xShelf  = xShelf(:);
        hShelf  = hShelf(:);
        bShelf  = bShelf(:);
        xSheet  = flipud(xSheet);
        hSheet  = flipud(hSheet);
        bSheet  = flipud(bSheet);
        hxSheet = flipud(hxSheet);
        HSheet = hSheet - bSheet;
        HShelf = hShelf - bShelf;
        xN = xShelf(end);
        if scaleXByGroundingLine
            xScale = xG;
        else
            xScale = 1;
        end
        result.modelKey = 'shear-dominated';
        result.beta = beta;
        result.epsilon = epsilon;
        result.Accumulation = Accumulation;
        result.Lambda = Lambda;
        result.n = n;
        result.Dn = 2^(1-n)/(n+2);
        result.xG = xG;
        result.HG = HG;
        result.hG = hG;
        result.bG = bG;
        result.xN = xN;
        result.xSheet = xSheet;
        result.HSheet = HSheet;
        result.hSheet = hSheet;
        result.bSheet = bSheet;
        result.hxSheet = hxSheet;
        result.xShelf = xShelf;
        result.HShelf = HShelf;
        result.hShelf = hShelf;
        result.bShelf = bShelf;
        result.xSheetPlot = xSheet./xScale;
        result.xShelfPlot = xShelf./xScale;
        result.xGPlot = xG./xScale;
        result.xNPlot = xN./xScale;
        result.xMaxPlot = xN./xScale;
    end

    function result = callExtensionalSolver( ...
        beta,epsilon,Accumulation,Lambda,n, shelfExtent,scaleXByGroundingLine)
        [xSheet,xShelf,hSheet,hShelf,bSheet,bShelf, xG,HG,hG,bG] = ...
            Extensionalicesheet( Lambda, beta, Accumulation, n, epsilon, ...
                shelfExtent);
        xSheet = xSheet(:);
        hSheet = hSheet(:);
        bSheet = bSheet(:);
        xShelf = xShelf(:);
        hShelf = hShelf(:);
        bShelf = bShelf(:);
        xSheet = flipud(xSheet);
        hSheet = flipud(hSheet);
        bSheet = flipud(bSheet);
        HSheet = hSheet - bSheet;
        HShelf = hShelf - bShelf;
        xN = xShelf(end);
        if scaleXByGroundingLine
            xScale = xG;
        else
            xScale = 1;
        end
        result.modelKey = 'extensional';
        result.beta = beta;
        result.epsilon = epsilon;
        result.Accumulation = Accumulation;
        result.Lambda = Lambda;
        result.n = n;
        result.Dn = 2^(1-n)/(n+2);
        result.xG = xG;
        result.HG = HG;
        result.hG = hG;
        result.bG = bG;
        result.xN = xN;
        result.xSheet = xSheet;
        result.HSheet = HSheet;
        result.hSheet = hSheet;
        result.bSheet = bSheet;
        result.xShelf = xShelf;
        result.HShelf = HShelf;
        result.hShelf = hShelf;
        result.bShelf = bShelf;
        result.xSheetPlot = xSheet./xScale;
        result.xShelfPlot = xShelf./xScale;
        result.xGPlot = xG./xScale;
        result.xNPlot = xN./xScale;
        result.xMaxPlot = xN./xScale;
    end

    function plotConstantDepthProfilesGrid( results,colours,arrowSettings, ...
        epsilon,nValues,lambdaValues, ...
        scaleXByGroundingLine,exportSettings,outputFilename)
        [nRows,nColumns] = size(results);
        yMinimum = inf;
        yMaximum = -inf;
        for row = 1:nRows
            for column = 1:nColumns
                result = results{row,column};
                yMinimum = min([ yMinimum;
                    result.bSheet(:);
                    result.bShelf(:)]);
                yMaximum = max([ yMaximum;
                    result.hSheet(:);
                    result.hShelf(:)]);
            end
        end
        ySpan = yMaximum - yMinimum;
        if ySpan <= 0
            ySpan = 1;
        end
        yBottom = yMinimum - 0.18*ySpan;
        yTop = yMaximum + 0.10*ySpan;
        maximumVelocityByColumn = zeros(1,nColumns);
        for column = 1:nColumns
            columnMaximum = 0;
            for row = 1:nRows
                columnMaximum = max( columnMaximum, ...
                    profileMaximumVelocity(results{row,column}));
            end
            maximumVelocityByColumn(column) = columnMaximum;
        end
        figureWidthCm = exportSettings.figureWidthCm;
        figureHeightCm = exportSettings.figureHeightCm;
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[2,2,figureWidthCm,figureHeightCm]);
        layout = tiledlayout( nRows,nColumns, 'TileSpacing','compact', ...
            'Padding','compact');
        layout.OuterPosition = [0.115,0.17,0.79,0.72];
        legendHandles = [];
        for row = 1:nRows
            for column = 1:nColumns
                result = results{row,column};
                ax = nexttile(layout);
                handles = plotIceProfile( ax, result, yBottom, yTop, ...
                    colours, maximumVelocityByColumn(column), arrowSettings);
                if row == nRows && column == 1
                    if scaleXByGroundingLine
                        xlabel(ax,'$x/x_G$','Interpreter','latex','FontSize',9);
                    else
                        xlabel(ax,'$x$','Interpreter','latex','FontSize',9);
                    end
                    ylabel(ax,'$z$','Interpreter','latex','FontSize',9);
                    ax.XLabel.Units = 'normalized';
                    ax.XLabel.Position = [0.5,-0.105,0];
                    ax.YLabel.Units = 'normalized';
                    ax.YLabel.Position = [-0.085,0.5,0];
                    ax.XTickLabelMode = 'auto';
                    ax.YTickLabelMode = 'auto';
                    ax.YTick = chooseIntegerTicks(yBottom,yTop,5);
                else
                    ax.XTickLabel = [];
                    ax.YTickLabel = [];
                end
                if row == 1 && column == 1
                    legendHandles = [ handles.ice, handles.bedrock, ...
                        handles.ocean, handles.velocity];
                end
            end
        end
        sharedLegend = legend( legendHandles, ...
            {'Ice','Bedrock','Ocean','Velocity'}, 'Interpreter','none', ...
            'FontSize',9, 'Orientation','horizontal', 'Box','off');
        sharedLegend.NumColumns = 4;
        sharedLegend.Layout.Tile = 'south';
        addOuterParameterAxes( figureHandle, nValues, lambdaValues);
        exportgraphics( figureHandle, resultsPath(outputFilename), 'Resolution',300);
    end

    function addOuterParameterAxes(figureHandle,nValues,lambdaValues)
        nLabels = arrayfun( @(n) num2str(n,'%.4g'), nValues, ...
            'UniformOutput',false);
        lambdaLabels = arrayfun( @(lambda) num2str(lambda,'%.4g'), ...
            lambdaValues, 'UniformOutput',false);
        parameterAxes = axes( 'Parent',figureHandle, ...
            'Position',[0.105,0.175,0.80,0.71], 'Color','none', ...
            'Box','off', 'XAxisLocation','bottom', 'YAxisLocation','left', ...
            'XLim',[0.5,numel(nValues)+0.5], ...
            'YLim',[0.5,numel(lambdaValues)+0.5], ...
            'XTick',1:numel(nValues), 'YTick',1:numel(lambdaValues), ...
            'XTickLabel',nLabels, 'YTickLabel',lambdaLabels, ...
            'TickDir','out', 'TickLength',[0.012,0.012], 'LineWidth',1.10, ...
            'FontSize',9, 'Layer','top', 'HitTest','off', ...
            'PickableParts','none');
        parameterAxes.TickLabelInterpreter = 'none';
        xlabel(parameterAxes,'$n$', ...
            'Interpreter','latex','FontWeight','bold','FontSize',9);
        ylabel(parameterAxes,'$\Lambda$', ...
            'Interpreter','latex','FontWeight','bold','FontSize',9);
        parameterAxes.XLabel.Units = 'normalized';
        parameterAxes.XLabel.Position = [0.5,-0.12,0];
        parameterAxes.YLabel.Units = 'normalized';
        parameterAxes.YLabel.Position = [-0.065,0.5,0];
        parameterAxesPosition = parameterAxes.Position;
        annotation(figureHandle,'arrow', ...
            [parameterAxesPosition(1)+parameterAxesPosition(3), ...
             parameterAxesPosition(1)+parameterAxesPosition(3)+0.02], ...
            [parameterAxesPosition(2), parameterAxesPosition(2)], ...
            'LineWidth',1.10,'HeadLength',6,'HeadWidth',6);
        annotation(figureHandle,'arrow', ...
            [parameterAxesPosition(1), parameterAxesPosition(1)], ...
            [parameterAxesPosition(2)+parameterAxesPosition(4), ...
             parameterAxesPosition(2)+parameterAxesPosition(4)+0.02], ...
            'LineWidth',1.10,'HeadLength',6,'HeadWidth',6);
    end

    function yTicks = chooseIntegerTicks(yMinimum,yMaximum,maxTickCount)
        if nargin < 3
            maxTickCount = 5;
        end
        startTick = ceil(yMinimum);
        endTick = floor(yMaximum);
        if startTick > endTick
            yTicks = round(linspace(yMinimum,yMaximum,min(3,maxTickCount)));
            yTicks = unique(yTicks);
            return;
        end
        integerTicks = startTick:endTick;
        if numel(integerTicks) <= maxTickCount
            yTicks = integerTicks;
            return;
        end
        step = ceil((numel(integerTicks)-1)/(maxTickCount-1));
        yTicks = integerTicks(1:step:end);
        if yTicks(end) ~= integerTicks(end)
            yTicks(end+1) = integerTicks(end);
        end
    end

    function handles = plotIceProfile( ax,result,yBottom,yTop,colours, ...
        maximumVelocity,arrowSettings)
        hold(ax,'on');
        xSheet = result.xSheetPlot(:);
        xShelf = result.xShelfPlot(:);
        hSheet = result.hSheet(:);
        bSheet = result.bSheet(:);
        hShelf = result.hShelf(:);
        bShelf = result.bShelf(:);
        bOcean = -result.beta.*result.xShelf(:);
        handles.bedrock = patch( ax, [xSheet;flipud(xSheet)], ...
            [bSheet;yBottom.*ones(size(xSheet))], colours.bedrock, ...
            'EdgeColor','none');
        patch( ax, [xShelf;flipud(xShelf)], ...
            [bOcean;yBottom.*ones(size(xShelf))], colours.bedrock, ...
            'EdgeColor','none');
        handles.ocean = patch( ax, [xShelf;flipud(xShelf)], ...
            [bShelf;flipud(bOcean)], colours.ocean, 'EdgeColor','none');
        handles.ice = patch( ax, [xSheet;flipud(xSheet)], ...
            [hSheet;flipud(bSheet)], colours.ice, 'EdgeColor','none');
        patch( ax, [xShelf;flipud(xShelf)], [hShelf;flipud(bShelf)], ...
            colours.ice, 'EdgeColor','none');
        handles.velocity = addVelocityProfiles( ax, result, ...
            maximumVelocity, arrowSettings);
        plot( ax,xSheet,hSheet, 'Color',colours.surface, 'LineWidth',1.5);
        plot( ax,xShelf,hShelf, 'Color',colours.surface, 'LineWidth',1.5);
        plot( ax,xSheet,bSheet, 'Color',colours.bed, 'LineWidth',1.5);
        plot( ax,xShelf,bShelf, 'Color',colours.surface, 'LineWidth',1.5);
        plot( ax,xShelf,bOcean, 'Color',colours.bed, 'LineWidth',1.5);
        xline( ax, result.xGPlot, '--k', 'LineWidth',0.9);
        yline( ax, 0, '--', 'Color',colours.sea, 'LineWidth',0.9);
        xlim(ax,[0,result.xMaxPlot]);
        ylim(ax,[yBottom,yTop]);
        set( ax, 'TickLabelInterpreter','latex', 'FontSize',9);
        box(ax,'on');
        hold(ax,'off');
    end

    function maximumVelocity = profileMaximumVelocity(result)
        etaSurface = ones(size(result.HSheet));
        if strcmp(result.modelKey,'shear-dominated')
            uSheet = velocityprofileSheardominatedsheet( etaSurface, ...
                result.HSheet, result.hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            uSheet = velocityprofileExtensionalsheet( etaSurface, ...
                result.xSheet, result.HSheet, result.Accumulation);
        else
            error('Unknown model key: %s.',result.modelKey);
        end
        etaShelf = ones(size(result.HShelf));
        uShelf = velocityprofileshelf( etaShelf, result.xShelf, ...
            result.HShelf, result.Accumulation);
        maximumVelocity = max([ uSheet(:);
            uShelf(:)]);
    end

    function arrowHandle = addVelocityProfiles( ...
        ax,result,maximumVelocity,settings)
        if maximumVelocity <= 0 || ~isfinite(maximumVelocity)
            error('maximumVelocity must be positive and finite.');
        end
        etaValues = linspace(0.08,0.92,settings.nVerticalLevels).';
        xSheetPosition = settings.sheetArrowXFraction * result.xG;
        if result.xGPlot ~= 0
            xScale = result.xG/result.xGPlot;
        else
            xScale = 1;
        end
        xSheetPlotPosition = xSheetPosition/xScale;
        HSheet = interp1( result.xSheet, result.HSheet, xSheetPosition, ...
            'pchip');
        bSheet = interp1( result.xSheet, result.bSheet, xSheetPosition, ...
            'pchip');
        zSheet = bSheet + etaValues.*HSheet;
        if strcmp(result.modelKey,'shear-dominated')
            hxSheet = interp1( result.xSheet, result.hxSheet, ...
                xSheetPosition, 'pchip');
            uSheet = velocityprofileSheardominatedsheet( etaValues, ...
                HSheet, hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            uSheet = velocityprofileExtensionalsheet( etaValues, ...
                xSheetPosition, HSheet, result.Accumulation);
        else
            error('Unknown model key: %s.',result.modelKey);
        end
        xShelfPosition = result.xG + settings.shelfArrowXFraction ...
              * (result.xN - result.xG);
        xShelfPlotPosition = xShelfPosition/xScale;
        HShelf = interp1( result.xShelf, result.HShelf, xShelfPosition, ...
            'pchip');
        bShelf = interp1( result.xShelf, result.bShelf, xShelfPosition, ...
            'pchip');
        zShelf = bShelf + etaValues.*HShelf;
        uShelf = velocityprofileshelf( etaValues, xShelfPosition, HShelf, ...
            result.Accumulation);
        arrowX = [ xSheetPlotPosition.*ones(size(etaValues));
            xShelfPlotPosition.*ones(size(etaValues))];
        arrowZ = [ zSheet;
            zShelf];
        velocityValues = [ uSheet(:);
            uShelf(:)];
        maximumArrowLength = settings.maxLengthFraction * result.xMaxPlot;
        plottedArrowLengths = maximumArrowLength .* velocityValues ...
            ./ maximumVelocity;
        arrowHandle = quiver( ax, arrowX, arrowZ, plottedArrowLengths, ...
            zeros(size(plottedArrowLengths)), 0, 'Color',settings.colour, ...
            'LineWidth',settings.lineWidth, ...
            'MaxHeadSize',settings.headSize, 'AutoScale','off');
    end

    function exportConstantSheetDepthParameters(results,csvFilename,texFilename)
    %EXPORTCONSTANTSHEETDEPTHPARAMETERS
    % Export the bed slope beta used for each constant-sheet-depth profile.
    % The LaTeX output is arranged to match the plotted 2-by-2 grid:
    %       columns -> rheological exponent n
    %       rows    -> ice-sheet model
    % The n-values and Lambda-values are read directly from results, so the
    % table automatically follows any later changes to the plotted cases.
        nRows = size(results,1);
        nColumns = size(results,2);
        if nRows ~= 2
            error('The compact LaTeX table expects exactly two model rows.');
        end
        % -------------------------------------------
        % CSV export: one row per model, one column per n-value
        % -------------------------------------------
        Model = strings(nRows,1);
        Lambda = zeros(nRows,1);
        bedSlopeMatrix = zeros(nRows,nColumns);
        for row = 1:nRows
            referenceResult = results{row,1};
            if strcmp(referenceResult.modelKey,'extensional')
                Model(row) = "Extensional ice sheet";
            elseif strcmp(referenceResult.modelKey,'shear-dominated')
                Model(row) = "Shear-dominated ice sheet";
            else
                Model(row) = string(referenceResult.modelKey);
            end
            Lambda(row) = referenceResult.Lambda;
            for column = 1:nColumns
                bedSlopeMatrix(row,column) = results{row,column}.beta;
            end
        end
        csvTable = table(Model,Lambda);
        for column = 1:nColumns
            nValue = results{1,column}.n;
            variableName = matlab.lang.makeValidName(sprintf('n_%g',nValue));
            csvTable.(variableName) = bedSlopeMatrix(:,column);
        end
        writetable(csvTable,resultsPath(csvFilename));
        % -------------------------------------------
        % Compact LaTeX table
        % -------------------------------------------
        fid = fopen(resultsPath(texFilename),'w');
        if fid == -1
            error('Could not open %s for writing.',texFilename);
        end
        cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid, ...
            '%% Generated automatically by constant_ice_sheet_depth_profiles.m\n');        fprintf(fid,'\\begin{tabular}{l');
        for column = 1:nColumns
            fprintf(fid,'c');
        end
        fprintf(fid,'}\n');
        fprintf(fid,'\\toprule\n');
        % One compact header row only.
        fprintf(fid,' ');
        for column = 1:nColumns
            nValue = results{1,column}.n;
            fprintf(fid,' & $n=%g$',nValue);
        end
        fprintf(fid,' \\\\\n');
        fprintf(fid,'\\midrule\n');
        for row = 1:nRows
            result = results{row,1};
            if strcmp(result.modelKey,'extensional')
                modelLabel = sprintf( ...
                    'Extensional ice sheet ($\\Lambda=%g$)', result.Lambda);
            elseif strcmp(result.modelKey,'shear-dominated')
                modelLabel = sprintf( ...
                    'Shear-dominated ice sheet ($\\Lambda=%g$)', ...
                    result.Lambda);
            else
                modelLabel = char(string(result.modelKey));
            end
            fprintf(fid,'%s',modelLabel);
            for column = 1:nColumns
                fprintf( fid, ' & %s', latexNumber(results{row,column}.beta));
            end
            fprintf(fid,' \\\\\n');
        end
        fprintf(fid,'\\bottomrule\n');
        fprintf(fid,'\\end{tabular}\n');
    end

    function textValue = latexNumber(value)
    %LATEXNUMBER Format a MATLAB number cleanly for a LaTeX table.
        if value == 0
            textValue = '$0$';
            return;
        end
        exponent = floor(log10(abs(value)));
        if exponent <= -3 || exponent >= 4
            coefficient = value/10^exponent;
            if abs(coefficient - 1) < 1e-12
                textValue = sprintf('$10^{%d}$',exponent);
            elseif abs(coefficient + 1) < 1e-12
                textValue = sprintf('$-10^{%d}$',exponent);
            else
                textValue = sprintf('$%.4g\\times10^{%d}$',coefficient,exponent);
            end
        else
            textValue = sprintf('$%.4g$',value);
        end
    end
end

function runConstantIceShelfDepthStudy()
    %% Constant ice-shelf-depth profiles
    % Compare shear-dominated and extensional solutions for which the
    % floating shelf has constant depth, for n = 1 and n = 3.
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    %% MODEL PARAMETERS
    epsilon = 0.1;
    LambdaShearDominated = 0.3;
    LambdaExtensional = 2;
    beta = 5;
    shelfExtent = 3;
    scaleXByGroundingLine = true;

    %% COLOURS
    colours.ice      = [0.78, 0.91, 1.00];
    colours.bedrock  = [0.48, 0.30, 0.18];
    colours.ocean    = [0.25, 0.62, 0.85];
    colours.surface  = [0.05, 0.25, 0.65];
    colours.bed      = [0.25, 0.25, 0.25];
    colours.sea      = [0.10, 0.35, 0.55];
    colours.velocity = [0.05, 0.05, 0.05];

    %% VELOCITY ARROW SETTINGS
    velocityArrowSettings.nVerticalLevels = 3;
    velocityArrowSettings.sheetArrowXFraction = 0.45;
    velocityArrowSettings.shelfArrowXFraction = 0.45;
    velocityArrowSettings.maxLengthFraction = 0.11;
    velocityArrowSettings.colour = colours.velocity;
    velocityArrowSettings.lineWidth = 0.85;
    velocityArrowSettings.headSize = 0.75;

    %% EXPORT SETTINGS
    exportSettings.figureWidthCm = 16;
    exportSettings.figureHeightCm = 8.0;

    %% EXTENSIONAL ICE SHEET: n = 1
    extensional1Result = callExtensionalConstantShelfSolver( beta, ...
        epsilon, LambdaExtensional, 1, shelfExtent, scaleXByGroundingLine);

    %% EXTENSIONAL ICE SHEET: n = 3
    extensional3Result = callExtensionalConstantShelfSolver( beta, ...
        epsilon, LambdaExtensional, 3, shelfExtent, scaleXByGroundingLine);

    %% SHEAR-DOMINATED ICE SHEET: n = 1
    shear1Result = callShearDominatedConstantShelfSolver( beta, epsilon, ...
        LambdaShearDominated, 1, shelfExtent, scaleXByGroundingLine);

    %% SHEAR-DOMINATED ICE SHEET: n = 3
    shear3Result = callShearDominatedConstantShelfSolver( beta, epsilon, ...
        LambdaShearDominated, 3, shelfExtent, scaleXByGroundingLine);

    %% ARRANGE RESULTS
    %  Top row:
    %       Extensional ice sheet
    %  Bottom row:
    %       Shear-dominated ice sheet
    %  Left column:
    %       n = 1
    %  Right column:
    %       n = 3
    results = { extensional1Result,extensional3Result; ...
        shear1Result,shear3Result};

    %% GENERATE PROFILE FIGURE
    plotConstantShelfDepthProfiles( results, colours, ...
        velocityArrowSettings, [1,3], ...
        [LambdaShearDominated,LambdaExtensional], scaleXByGroundingLine, ...
        exportSettings, 'constant_ice_shelf_depth_profiles.png');

    %% CONSTANT-SHELF PROFILE HELPERS

    function result = callShearDominatedConstantShelfSolver( ...
        beta,epsilon,Lambda,n, shelfExtent,scaleXByGroundingLine)
        [ xSheet,xShelf,hSheet,hShelf, bSheet,bShelf,hxSheet,xG, ...
            HG,hG,bG,Accumulation] = SheardominatedicesheetConstantShelf( ...
                Lambda, beta, n, epsilon, shelfExtent);
        result = packageResult( 'shear-dominated', ...
            ['$\mathrm{Shear\!-\!dominated},\ n=', num2str(n,'%.4g'),'$'], ...
            beta, epsilon, Accumulation, Lambda, n, xSheet, xShelf, ...
            hSheet, hShelf, bSheet, bShelf, hxSheet, xG, HG, hG, bG, ...
            scaleXByGroundingLine);
    end

    %% EXTENSIONAL CONSTANT-SHELF WRAPPER

    function result = callExtensionalConstantShelfSolver( ...
        beta,epsilon,Lambda,n, shelfExtent,scaleXByGroundingLine)
        [ xSheet,xShelf,hSheet,hShelf, bSheet,bShelf,xG, ...
            HG,hG,bG,Accumulation] = ExtensionalicesheetConstantShelf( ...
                Lambda, beta, n, epsilon, shelfExtent);
        hxSheet = [];
        result = packageResult( 'extensional', ...
            ['$\mathrm{Extensional},\ n=', num2str(n,'%.4g'),'$'], beta, ...
            epsilon, Accumulation, Lambda, n, xSheet, xShelf, hSheet, ...
            hShelf, bSheet, bShelf, hxSheet, xG, HG, hG, bG, ...
            scaleXByGroundingLine);
    end

    %% PACKAGE RESULT

    function result = packageResult( modelKey,titleLabel, ...
        beta,epsilon,Accumulation,Lambda,n, ...
        xSheet,xShelf,hSheet,hShelf,bSheet,bShelf,hxSheet, ...
        xG,HG,hG,bG,scaleXByGroundingLine)
        % Store profile outputs as column vectors.
        xSheet = xSheet(:);
        hSheet = hSheet(:);
        bSheet = bSheet(:);
        xShelf = xShelf(:);
        hShelf = hShelf(:);
        bShelf = bShelf(:);
        % The sheet is integrated backwards from x_G to x = 0.
        % Reverse it so that x increases from left to right.
        xSheet = flipud(xSheet);
        hSheet = flipud(hSheet);
        bSheet = flipud(bSheet);
        if ~isempty(hxSheet)
            hxSheet = flipud(hxSheet(:));
        end
        % Ice thickness.
        HSheet = hSheet - bSheet;
        HShelf = hShelf - bShelf;
        % Shelf end position.
        xN = xShelf(end);
        % Horizontal plotting scale.
        if scaleXByGroundingLine
            xScale = xG;
        else
            xScale = 1;
        end
        % Model information.
        result.modelKey = modelKey;
        result.titleLabel = titleLabel;
        result.beta = beta;
        result.epsilon = epsilon;
        result.Accumulation = Accumulation;
        result.Lambda = Lambda;
        result.n = n;
        % Grounding-line information.
        result.xG = xG;
        result.HG = HG;
        result.hG = hG;
        result.bG = bG;
        % Shelf end.
        result.xN = xN;
        % Sheet solution.
        result.xSheet = xSheet;
        result.HSheet = HSheet;
        result.hSheet = hSheet;
        result.bSheet = bSheet;
        result.hxSheet = hxSheet;
        % Shelf solution.
        result.xShelf = xShelf;
        result.HShelf = HShelf;
        result.hShelf = hShelf;
        result.bShelf = bShelf;
        % Coordinates used in the figure.
        result.xSheetPlot = xSheet./xScale;
        result.xShelfPlot = xShelf./xScale;
        result.xGPlot = xG./xScale;
        result.xNPlot = xN./xScale;
        result.xMaxPlot = xN./xScale;
    end

    %% PLOT CONSTANT-SHELF-DEPTH PROFILES

    function plotConstantShelfDepthProfiles( ...
        results,colours,arrowSettings, nValues,lambdaValues, ...
        scaleXByGroundingLine,exportSettings,outputFilename)
        [nRows,nColumns] = size(results);
        yMinimum = inf;
        yMaximum = -inf;
        maximumVelocity = 0;
        % -------------------------------------------
        % Find common vertical limits and velocity scaling
        % -------------------------------------------
        for k = 1:numel(results)
            result = results{k};
            yMinimum = min([ yMinimum;
                result.bSheet(:);
                result.bShelf(:)]);
            yMaximum = max([ yMaximum;
                result.hSheet(:);
                result.hShelf(:)]);
            maximumVelocity = max( maximumVelocity, ...
                profileMaximumVelocity(result));
        end
        ySpan = yMaximum - yMinimum;
        if ySpan <= 0
            ySpan = 1;
        end
        yBottom = yMinimum - 0.18*ySpan;
        yTop = yMaximum + 0.10*ySpan;
        % -------------------------------------------
        % Create figure
        % -------------------------------------------
        figureWidthCm = exportSettings.figureWidthCm;
        figureHeightCm = exportSettings.figureHeightCm;
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[ 2, 2, figureWidthCm, figureHeightCm]);
        % -------------------------------------------
        % Tiled layout
        % -------------------------------------------
        layout = tiledlayout( nRows, nColumns, 'TileSpacing','compact', ...
            'Padding','compact');
        layout.OuterPosition = [ 0.115, 0.17, 0.79, 0.72];
        legendHandles = [];
        % -------------------------------------------
        % Plot each panel
        % -------------------------------------------
        for row = 1:nRows
            for column = 1:nColumns
                result = results{row,column};
                ax = nexttile(layout);
                handles = plotIceProfile( ax, result, yBottom, yTop, ...
                    colours, maximumVelocity, arrowSettings);
                % -------------------------------------------
                % Only show numerical axis labels on lower-left panel
                % -------------------------------------------
                if row == nRows && column == 1
                    if scaleXByGroundingLine
                        xlabel( ax, '$x/x_G$', 'Interpreter','latex', ...
                            'FontSize',9);
                    else
                        xlabel( ax, '$x$', 'Interpreter','latex', ...
                            'FontSize',9);
                    end
                    ylabel( ax, '$z$', 'Interpreter','latex', 'FontSize',9);
                    ax.XLabel.Units = 'normalized';
                    ax.XLabel.Position = [0.5,-0.105,0];
                    ax.YLabel.Units = 'normalized';
                    ax.YLabel.Position = [-0.085,0.5,0];
                    ax.XTickLabelMode = 'auto';
                    ax.YTickLabelMode = 'auto';
                    ax.YTick = chooseIntegerTicks( yBottom, yTop, 5);
                else
                    ax.XTickLabel = [];
                    ax.YTickLabel = [];
                end
                % -------------------------------------------
                % Save handles for shared legend
                % -------------------------------------------
                if row == 1 && column == 1
                    legendHandles = [ handles.ice, handles.bedrock, ...
                        handles.ocean, handles.velocity];
                end
            end
        end
        % -------------------------------------------
        % Shared legend
        % -------------------------------------------
        sharedLegend = legend( legendHandles, { 'Ice', 'Bedrock', 'Ocean', ...
                'Velocity'}, 'Interpreter','none', 'FontSize',9, ...
            'Orientation','horizontal', 'Box','off');
        sharedLegend.NumColumns = 4;
        sharedLegend.Layout.Tile = 'south';
        % -------------------------------------------
        % Add outer n and Lambda axes
        % -------------------------------------------
        addOuterParameterAxes( figureHandle, nValues, lambdaValues);
        % -------------------------------------------
        % Export figure
        % -------------------------------------------
        exportgraphics( figureHandle, resultsPath(outputFilename), 'Resolution',300);
    end

    %% OUTER PARAMETER AXES

    function addOuterParameterAxes( figureHandle,nValues,lambdaValues)
        nLabels = arrayfun( @(n) num2str(n,'%.4g'), nValues, ...
            'UniformOutput',false);
        lambdaLabels = arrayfun( @(lambda) num2str(lambda,'%.4g'), ...
            lambdaValues, 'UniformOutput',false);
        parameterAxes = axes( 'Parent',figureHandle, ...
            'Position',[0.105,0.175,0.80,0.71], 'Color','none', ...
            'Box','off', 'XAxisLocation','bottom', 'YAxisLocation','left', ...
            'XLim',[0.5,numel(nValues)+0.5], ...
            'YLim',[0.5,numel(lambdaValues)+0.5], ...
            'XTick',1:numel(nValues), 'YTick',1:numel(lambdaValues), ...
            'XTickLabel',nLabels, 'YTickLabel',lambdaLabels, ...
            'TickDir','out', 'TickLength',[0.012,0.012], 'LineWidth',1.10, ...
            'FontSize',9, 'Layer','top', 'HitTest','off', ...
            'PickableParts','none');
        parameterAxes.TickLabelInterpreter = 'none';
        xlabel( parameterAxes, '$n$', 'Interpreter','latex', ...
            'FontWeight','bold', 'FontSize',9);
        ylabel( parameterAxes, '$\Lambda$', 'Interpreter','latex', ...
            'FontWeight','bold', 'FontSize',9);
        parameterAxes.XLabel.Units = 'normalized';
        parameterAxes.XLabel.Position = [0.5,-0.12,0];
        parameterAxes.YLabel.Units = 'normalized';
        parameterAxes.YLabel.Position = [-0.065,0.5,0];
        parameterAxesPosition = parameterAxes.Position;
        % -------------------------------------------
        % Horizontal axis arrow
        % -------------------------------------------
        annotation( figureHandle, 'arrow', [ parameterAxesPosition(1) ...
                + parameterAxesPosition(3), parameterAxesPosition(1) ...
                + parameterAxesPosition(3) + 0.02], [ ...
                parameterAxesPosition(2), parameterAxesPosition(2)], ...
            'LineWidth',1.10, 'HeadLength',6, 'HeadWidth',6);
        % -------------------------------------------
        % Vertical axis arrow
        % -------------------------------------------
        annotation( figureHandle, 'arrow', [ parameterAxesPosition(1), ...
                parameterAxesPosition(1)], [ parameterAxesPosition(2) ...
                + parameterAxesPosition(4), parameterAxesPosition(2) ...
                + parameterAxesPosition(4) + 0.02], 'LineWidth',1.10, ...
            'HeadLength',6, 'HeadWidth',6);
    end

    %% CHOOSE INTEGER Y TICKS

    function yTicks = chooseIntegerTicks( yMinimum,yMaximum,maxTickCount)
        if nargin < 3
            maxTickCount = 5;
        end
        startTick = ceil(yMinimum);
        endTick = floor(yMaximum);
        if startTick > endTick
            yTicks = round( linspace( yMinimum, yMaximum, ...
                    min(3,maxTickCount)));
            yTicks = unique(yTicks);
            return;
        end
        integerTicks = startTick:endTick;
        if numel(integerTicks) <= maxTickCount
            yTicks = integerTicks;
            return;
        end
        step = ceil( (numel(integerTicks)-1) /(maxTickCount-1));
        yTicks = integerTicks(1:step:end);
        if yTicks(end) ~= integerTicks(end)
            yTicks(end+1) = integerTicks(end);
        end
    end

    %% PLOT AN INDIVIDUAL ICE PROFILE

    function handles = plotIceProfile( ax,result,yBottom,yTop,colours, ...
        maximumVelocity,arrowSettings)
        hold(ax,'on');
        xSheet = result.xSheetPlot(:);
        xShelf = result.xShelfPlot(:);
        hSheet = result.hSheet(:);
        bSheet = result.bSheet(:);
        hShelf = result.hShelf(:);
        bShelf = result.bShelf(:);
        % Bedrock continues beneath the ocean according to
        % b = -beta*x.
        bOcean = -result.beta .* result.xShelf(:);
        % -------------------------------------------
        % Bedrock beneath grounded ice
        % -------------------------------------------
        handles.bedrock = patch( ax, [ xSheet; flipud(xSheet)], [ bSheet; ...
                yBottom.*ones(size(xSheet))], colours.bedrock, ...
            'EdgeColor','none');
        % -------------------------------------------
        % Bedrock beneath ocean
        % -------------------------------------------
        patch( ax, [ xShelf; flipud(xShelf)], [ bOcean; ...
                yBottom.*ones(size(xShelf))], colours.bedrock, ...
            'EdgeColor','none');
        % -------------------------------------------
        % Ocean
        % -------------------------------------------
        handles.ocean = patch( ax, [ xShelf; flipud(xShelf)], [ bShelf; ...
                flipud(bOcean)], colours.ocean, 'EdgeColor','none');
        % -------------------------------------------
        % Grounded ice
        % -------------------------------------------
        handles.ice = patch( ax, [ xSheet; flipud(xSheet)], [ hSheet; ...
                flipud(bSheet)], colours.ice, 'EdgeColor','none');
        % -------------------------------------------
        % Floating ice shelf
        % -------------------------------------------
        patch( ax, [ xShelf; flipud(xShelf)], [ hShelf; flipud(bShelf)], ...
            colours.ice, 'EdgeColor','none');
        % -------------------------------------------
        % Velocity arrows
        % -------------------------------------------
        handles.velocity = addVelocityProfiles( ax, result, ...
            maximumVelocity, arrowSettings);
        % -------------------------------------------
        % Ice surface
        % -------------------------------------------
        plot( ax, xSheet, hSheet, 'Color',colours.surface, 'LineWidth',1.5);
        plot( ax, xShelf, hShelf, 'Color',colours.surface, 'LineWidth',1.5);
        % -------------------------------------------
        % Grounded bed
        % -------------------------------------------
        plot( ax, xSheet, bSheet, 'Color',colours.bed, 'LineWidth',1.5);
        % -------------------------------------------
        % Lower ice-shelf surface
        % -------------------------------------------
        plot( ax, xShelf, bShelf, 'Color',colours.surface, 'LineWidth',1.5);
        % -------------------------------------------
        % Bed beneath ocean
        % -------------------------------------------
        plot( ax, xShelf, bOcean, 'Color',colours.bed, 'LineWidth',1.5);
        % -------------------------------------------
        % Grounding line
        % -------------------------------------------
        xline( ax, result.xGPlot, '--k', 'LineWidth',0.9);
        % -------------------------------------------
        % Sea level
        % -------------------------------------------
        yline( ax, 0, '--', 'Color',colours.sea, 'LineWidth',0.9);
        % -------------------------------------------
        % Axes
        % -------------------------------------------
        xlim( ax, [0,result.xMaxPlot]);
        ylim( ax, [yBottom,yTop]);
        set( ax, 'TickLabelInterpreter','latex', 'FontSize',9);
        box(ax,'on');
        hold(ax,'off');
    end

    %% MAXIMUM VELOCITY IN EACH PROFILE

    function maximumVelocity = profileMaximumVelocity(result)
        etaSurface = ones(size(result.HSheet));
        % -------------------------------------------
        % Ice-sheet velocity
        % -------------------------------------------
        if strcmp(result.modelKey,'shear-dominated')
            uSheet = velocityprofileSheardominatedsheet( etaSurface, ...
                result.HSheet, result.hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            uSheet = velocityprofileExtensionalsheet( etaSurface, ...
                result.xSheet, result.HSheet, result.Accumulation);
        else
            error( 'Unknown model key: %s.', result.modelKey);
        end
        % -------------------------------------------
        % Ice-shelf velocity
        % -------------------------------------------
        etaShelf = ones(size(result.HShelf));
        uShelf = velocityprofileshelf( etaShelf, result.xShelf, ...
            result.HShelf, result.Accumulation);
        % -------------------------------------------
        % Maximum velocity
        % -------------------------------------------
        maximumVelocity = max([ uSheet(:);
            uShelf(:)]);
    end

    %% ADD VELOCITY PROFILES

    function arrowHandle = addVelocityProfiles( ...
        ax,result,maximumVelocity,settings)
        if maximumVelocity <= 0 || ~isfinite(maximumVelocity)
            error( 'maximumVelocity must be positive and finite.');
        end
        % -------------------------------------------
        % Vertical positions of arrows
        % -------------------------------------------
        etaValues = linspace( 0.08, 0.92, settings.nVerticalLevels).';
        % -------------------------------------------
        % SHEET VELOCITY
        % -------------------------------------------
        xSheetPosition = settings.sheetArrowXFraction * result.xG;
        if result.xGPlot ~= 0
            xScale = result.xG / result.xGPlot;
        else
            xScale = 1;
        end
        xSheetPlotPosition = xSheetPosition / xScale;
        % -------------------------------------------
        % Interpolate thickness and bed elevation
        % -------------------------------------------
        HSheet = interp1( result.xSheet, result.HSheet, xSheetPosition, ...
            'pchip');
        bSheet = interp1( result.xSheet, result.bSheet, xSheetPosition, ...
            'pchip');
        zSheet = bSheet + etaValues.*HSheet;
        % -------------------------------------------
        % Calculate sheet velocity
        % -------------------------------------------
        if strcmp(result.modelKey,'shear-dominated')
            hxSheet = interp1( result.xSheet, result.hxSheet, ...
                xSheetPosition, 'pchip');
            uSheet = velocityprofileSheardominatedsheet( etaValues, ...
                HSheet, hxSheet, result.Lambda, result.n);
        elseif strcmp(result.modelKey,'extensional')
            uSheet = velocityprofileExtensionalsheet( etaValues, ...
                xSheetPosition, HSheet, result.Accumulation);
        else
            error( 'Unknown model key: %s.', result.modelKey);
        end
        % -------------------------------------------
        % SHELF VELOCITY
        % -------------------------------------------
        xShelfPosition = result.xG + settings.shelfArrowXFraction ...
            * (result.xN - result.xG);
        xShelfPlotPosition = xShelfPosition / xScale;
        % -------------------------------------------
        % Interpolate shelf geometry
        % -------------------------------------------
        HShelf = interp1( result.xShelf, result.HShelf, xShelfPosition, ...
            'pchip');
        bShelf = interp1( result.xShelf, result.bShelf, xShelfPosition, ...
            'pchip');
        zShelf = bShelf + etaValues.*HShelf;
        % -------------------------------------------
        % Shelf velocity
        % -------------------------------------------
        uShelf = velocityprofileshelf( etaValues, xShelfPosition, HShelf, ...
            result.Accumulation);
        % -------------------------------------------
        % DRAW VELOCITY ARROWS
        % -------------------------------------------
        arrowX = [ xSheetPlotPosition.*ones(size(etaValues));
            xShelfPlotPosition.*ones(size(etaValues))];
        arrowZ = [ zSheet;
            zShelf];
        velocityValues = [ uSheet(:);
            uShelf(:)];
        maximumArrowLength = settings.maxLengthFraction * result.xMaxPlot;
        plottedArrowLengths = maximumArrowLength .* velocityValues ...
            ./ maximumVelocity;
        arrowHandle = quiver( ax, arrowX, arrowZ, plottedArrowLengths, ...
            zeros(size(plottedArrowLengths)), 0, 'Color',settings.colour, ...
            'LineWidth',settings.lineWidth, ...
            'MaxHeadSize',settings.headSize, 'AutoScale','off');
    end
end

function runGroundingLineSweepStudy()
    %% GROUNDING-LINE PARAMETER SWEEPS
    %  This routine evaluates the grounding-line position x_G while varying
    %  Lambda, bed slope beta, accumulation rate and epsilon. The
    %  shear-dominated and extensional models are evaluated for the selected
    %  rheological exponents, and theoretical power-law slope triangles are
    %  added for comparison with the numerical trends.
    %  Separate plots are exported for each varied parameter and each value
    %  of n.

    %% Fixed reference parameters
    fixed.beta                = 1;
    fixed.accumulationRate = 1;
    fixed.epsilon          = 1e-6;
    fixed.Lambda           = 0.1;
    % Rheological exponents to plot.
    nPlotValues = [1, 3, 5];

    %% Parameter ranges
    % Lambda on the x-axis, with beta represented by colour.
    LambdaForLambdaPlot = logspace(-4, 6, 250);
    betaColour             = logspace(-1, 1, 7);
    % Lambda represented by colour in the other panels.
    LambdaColour      = logspace(-4, 6, 7);
    LambdaColourDense = logspace(-4, 6, 7);
    % x-axis values for the remaining panels.
    betaValues            = logspace(-4, 2, 50);
    accumulationValues = logspace(-4, 6, 60);
    epsilonValues      = logspace(-6, -1, 50);
    % Parameters used only for the epsilon sweep.
    betaEpsilon            = 1e-6;
    accumulation_epsilon = 1e-2;

    %% Extensional ice-sheet model plotting settings
    extensional.enabled   = true;
    extensional.lineWidth = 1.35;
    extensional.lineStyle = '--';

    %% Theoretical power-law triangle settings
    % Placement parameters control the position and orientation of slope triangles.
    % Decreasing Lambda sweeps.
    trianglePlacement.Lambda.shearDominated.xFraction = ...
        [0.68,0.80; 0.68,0.80; 0.68,0.80];
    trianglePlacement.Lambda.shearDominated.relativeTo = ...
        {'below','below','below'};
    trianglePlacement.Lambda.shearDominated.gapFraction  = [0.050,0.050,0.050];
    trianglePlacement.Lambda.shearDominated.orientation = ...
        {'standard','standard','standard'};
    trianglePlacement.Lambda.extensional.xFraction = ...
        [0.78,0.90; 0.78,0.90; 0.78,0.90];
    trianglePlacement.Lambda.extensional.relativeTo = ...
        {'above','above','above'};
    trianglePlacement.Lambda.extensional.gapFraction     = [0.050,0.050,0.050];
    trianglePlacement.Lambda.extensional.orientation = ...
        {'opposite','opposite','opposite'};
    % Decreasing beta sweeps.
    trianglePlacement.beta.shearDominated.xFraction = ...
        [0.66,0.78; 0.67,0.79; 0.67,0.79];
    trianglePlacement.beta.shearDominated.relativeTo = ...
        {'below','below','below'};
    trianglePlacement.beta.shearDominated.gapFraction = ...
        [0.050,0.050,0.050];
    trianglePlacement.beta.shearDominated.orientation = ...
        {'standard','standard','standard'};
    trianglePlacement.beta.extensional.xFraction = ...
        [0.78,0.90; 0.78,0.90; 0.78,0.90];
    trianglePlacement.beta.extensional.relativeTo = ...
        {'above','above','above'};
    trianglePlacement.beta.extensional.gapFraction = ...
        [0.050,0.050,0.050];
    trianglePlacement.beta.extensional.orientation = ...
        {'opposite','opposite','opposite'};
    % Increasing accumulation-rate sweeps.
    % These placements are selected from the plotted figures so that the
    % shear-dominated and extensional slope triangles remain distinct from
    % the corresponding curve families. Rows correspond to n = [1, 3, 5].
    trianglePlacement.accumulationRate.shearDominated.xFraction = ...
        [0.76,0.84; 0.76,0.84; 0.76,0.84];
    % These positions place the shear-dominated triangles below the
    % asymptotic solid-curve family without obscuring the plotted curves.
    trianglePlacement.accumulationRate.shearDominated.verticalAnchorFraction = ...
        [0.34,0.51,0.55];
    trianglePlacement.accumulationRate.shearDominated.relativeTo = ...
        {'below','below','below'};
    trianglePlacement.accumulationRate.shearDominated.gapFraction = ...
        [0.030,0.030,0.030];
    trianglePlacement.accumulationRate.shearDominated.orientation = ...
        {'opposite','opposite','opposite'};
    trianglePlacement.accumulationRate.shearDominated.sizeScale = ...
        [1.22,1.22,1.22];
    trianglePlacement.accumulationRate.extensional.xFraction = ...
        [0.81,0.89; 0.81,0.89; 0.81,0.89];
    % The extensional triangles are placed above the dashed family while
    % retaining sufficient space for the triangle labels.
    trianglePlacement.accumulationRate.extensional.verticalAnchorFraction = ...
        [0.87,0.86,0.86];
    trianglePlacement.accumulationRate.extensional.relativeTo = ...
        {'above','above','above'};
    trianglePlacement.accumulationRate.extensional.gapFraction = ...
        [0.045,0.045,0.045];
    trianglePlacement.accumulationRate.extensional.orientation = ...
        {'standard','standard','standard'};
    trianglePlacement.accumulationRate.extensional.sizeScale = ...
        [1.18,1.18,1.18];
    % Decreasing epsilon sweeps. Only the shear-dominated triangle is drawn.
    trianglePlacement.epsilon.shearDominated.xFraction = ...
        [0.68,0.80; 0.68,0.80; 0.68,0.80];
    trianglePlacement.epsilon.shearDominated.relativeTo = ...
        {'below','below','below'};
    trianglePlacement.epsilon.shearDominated.gapFraction  = [0.045,0.045,0.045];
    trianglePlacement.epsilon.shearDominated.orientation = ...
        {'standard','standard','standard'};
    trianglePlacement.epsilon.extensional.xFraction = ...
        [0.82,0.94; 0.82,0.94; 0.82,0.94];
    trianglePlacement.epsilon.extensional.relativeTo = ...
        {'above','above','above'};
    trianglePlacement.epsilon.extensional.gapFraction     = [0.050,0.050,0.050];
    trianglePlacement.epsilon.extensional.orientation = ...
        {'opposite','opposite','opposite'};
    triangleStyle.shearDominated.Color = [0.0, 0.0, 0.0];
    triangleStyle.shearDominated.LineStyle = '-';
    triangleStyle.shearDominated.LineWidth = 1.15;
    triangleStyle.extensional.Color = [0.65, 0.10, 0.10];
    triangleStyle.extensional.LineStyle = '--';
    triangleStyle.extensional.LineWidth = 1.15;
    % Extensional-model power-law exponents used for the slope triangles.
    extensionalExponents.Lambda           = [-1/5, -1/19, -1/41];
    extensionalExponents.beta                = [-1, -1, -1];
    extensionalExponents.accumulationRate = [2/3, 4/15, 6/35];
    extensionalExponents.epsilon          = [-1, -1, -1];

    %% Figure export settings
    exportStyle.figureWidthCm       = 9.2;
    exportStyle.figureHeightCm      = 7.2;
    exportStyle.axesFontSizePt      = 9;
    exportStyle.labelFontSizePt     = 9;
    exportStyle.colourBarFontSizePt = 9;
    exportStyle.triangleFontSizePt   = 9;
    exportStyle.lineWidth            = 1.0;
    exportStyle.axesPosition       = [0.12, 0.18, 0.66, 0.74];
    exportStyle.colourBarPosition   = [0.82, 0.22, 0.025, 0.64];
    exportStyle.lowerYLimitPadding = 0.08;
    exportStyle.upperYLimitPadding = 0.06;
    exportStyle.accumulationUpperYLimitPadding = 0.12;
    exportResolution = 300;

    %% GENERATE PARAMETER SWEEPS FOR EACH n
    for nIndex = 1:numel(nPlotValues)
        nCurrent = nPlotValues(nIndex);
        fixedN = fixed;
        fixedN.n = nCurrent;
        % -------------------------------------------
        % a. Variation with Lambda, colour = beta
        % -------------------------------------------
        xGShearDominated = groundingLineSweep( 'shearDominated', ...
            'Lambda',LambdaForLambdaPlot, 'beta',betaColour, fixedN);
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[2,2, ...
            exportStyle.figureWidthCm,exportStyle.figureHeightCm]);
        ax = axes(figureHandle);
        loglog(ax,LambdaForLambdaPlot,xGShearDominated);
        formatSweepAxes(ax,'$\Lambda$',exportStyle);
        colourGradient = applyColourGradient(ax,betaColour,'$\beta$',exportStyle);
        lambdaExponent = expectedPowerLawExponent('Lambda',nCurrent);
        xGExtensional = [];
        if extensional.enabled
            xGExtensional = groundingLineSweep( 'extensional', ...
                'Lambda',LambdaForLambdaPlot, 'beta',betaColour, fixedN);
            plotExtensionalFamily( ax,LambdaForLambdaPlot,xGExtensional, ...
                colourGradient,extensional.lineWidth,extensional.lineStyle);
        end
        setSweepAxisLimits( ...
            ax,LambdaForLambdaPlot,xGShearDominated, ...
            xGExtensional,exportStyle);
        shearPlacement = selectTrianglePlacement( ...
            trianglePlacement,'Lambda','shearDominated',nCurrent,nPlotValues);
        addPowerLawTriangle( ...
            ax,LambdaForLambdaPlot,xGShearDominated,lambdaExponent, ...
            shearPlacement, triangleStyle.shearDominated,exportStyle);
        if extensional.enabled
            extensionalExponent = configuredExtensionalExponent( ...
                extensionalExponents,'Lambda',nCurrent,nPlotValues);
            if isfinite(extensionalExponent)
                extensionalPlacement = selectTrianglePlacement( ...
                    trianglePlacement,'Lambda','extensional', ...
                    nCurrent,nPlotValues);
                addPowerLawTriangle( ...
                    ax,LambdaForLambdaPlot,xGExtensional,extensionalExponent, ...
                    extensionalPlacement, ...
                    triangleStyle.extensional,exportStyle);
            end
        end
        outputFilename = sprintf( 'grounding_line_lambda_sweep_n%g.png',nCurrent);
        exportgraphics( figureHandle,resultsPath(outputFilename), ...
            'Resolution',exportResolution);
        % -------------------------------------------
        % b. Variation with beta, colour = Lambda
        % -------------------------------------------
        xGShearDominated = groundingLineSweep( 'shearDominated', ...
            'beta',betaValues, 'Lambda',LambdaColourDense, fixedN);
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[2,2, ...
            exportStyle.figureWidthCm,exportStyle.figureHeightCm]);
        ax = axes(figureHandle);
        loglog(ax,betaValues,xGShearDominated);
        formatSweepAxes(ax,'$\beta$',exportStyle);
        colourGradient = applyColourGradient( ...
            ax,LambdaColourDense,'$\Lambda$',exportStyle);
        betaExponent = expectedPowerLawExponent('beta',nCurrent);
        xGExtensional = [];
        if extensional.enabled
            xGExtensional = groundingLineSweep( 'extensional', ...
                'beta',betaValues, 'Lambda',LambdaColourDense, fixedN);
            plotExtensionalFamily( ax,betaValues,xGExtensional, ...
                colourGradient,extensional.lineWidth,extensional.lineStyle);
        end
        setSweepAxisLimits( ...
            ax,betaValues,xGShearDominated, ...
            xGExtensional,exportStyle);
        shearPlacement = selectTrianglePlacement( ...
            trianglePlacement,'beta','shearDominated',nCurrent,nPlotValues);
        addPowerLawTriangle( ax,betaValues,xGShearDominated,betaExponent, ...
            shearPlacement, triangleStyle.shearDominated,exportStyle);
        if extensional.enabled
            extensionalExponent = configuredExtensionalExponent( ...
                extensionalExponents,'beta',nCurrent,nPlotValues);
            if isfinite(extensionalExponent)
                extensionalPlacement = selectTrianglePlacement( ...
                    trianglePlacement,'beta','extensional',nCurrent,nPlotValues);
                addPowerLawTriangle( ...
                    ax,betaValues,xGExtensional,extensionalExponent, ...
                    extensionalPlacement, ...
                    triangleStyle.extensional,exportStyle);
            end
        end
        outputFilename = sprintf( 'grounding_line_beta_sweep_n%g.png',nCurrent);
        exportgraphics( figureHandle,resultsPath(outputFilename), ...
            'Resolution',exportResolution);
        % -------------------------------------------
        % c. Variation with accumulation, colour = Lambda
        % -------------------------------------------
        xGShearDominated = groundingLineSweep( 'shearDominated', ...
            'accumulationRate',accumulationValues, 'Lambda',LambdaColour, ...
            fixedN);
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[2,2, ...
            exportStyle.figureWidthCm,exportStyle.figureHeightCm]);
        ax = axes(figureHandle);
        loglog(ax,accumulationValues,xGShearDominated);
        formatSweepAxes(ax,'$\mathcal{A}$',exportStyle);
        colourGradient = applyColourGradient( ...
            ax,LambdaColour,'$\Lambda$',exportStyle);
        accumulationExponent = expectedPowerLawExponent( ...
            'accumulationRate',nCurrent);
        xGExtensional = [];
        if extensional.enabled
            xGExtensional = groundingLineSweep( 'extensional', ...
                'accumulationRate',accumulationValues, ...
                'Lambda',LambdaColour, fixedN);
            plotExtensionalFamily( ax,accumulationValues,xGExtensional, ...
                colourGradient,extensional.lineWidth,extensional.lineStyle);
        end
        accumulationExportStyle = exportStyle;
        accumulationExportStyle.upperYLimitPadding = ...
            exportStyle.accumulationUpperYLimitPadding;
        setSweepAxisLimits( ...
            ax,accumulationValues,xGShearDominated, ...
            xGExtensional,accumulationExportStyle);
        shearPlacement = selectTrianglePlacement( ...
            trianglePlacement,'accumulationRate','shearDominated', ...
            nCurrent,nPlotValues);
        addPowerLawTriangle( ...
            ax,accumulationValues,xGShearDominated,accumulationExponent, ...
            shearPlacement, triangleStyle.shearDominated,exportStyle);
        if extensional.enabled
            extensionalExponent = configuredExtensionalExponent( ...
                extensionalExponents,'accumulationRate',nCurrent,nPlotValues);
            if isfinite(extensionalExponent)
                extensionalPlacement = selectTrianglePlacement( ...
                    trianglePlacement,'accumulationRate','extensional', ...
                    nCurrent,nPlotValues);
                addPowerLawTriangle( ...
                    ax,accumulationValues,xGExtensional,extensionalExponent, ...
                    extensionalPlacement, ...
                    triangleStyle.extensional,exportStyle);
            end
        end
        outputFilename = sprintf( 'grounding_line_accumulation_sweep_n%g.png',nCurrent);
        exportgraphics( figureHandle,resultsPath(outputFilename), ...
            'Resolution',exportResolution);
        % -------------------------------------------
        % d. Variation with epsilon, colour = Lambda
        % -------------------------------------------
        fixedEpsilon = fixedN;
        fixedEpsilon.beta = betaEpsilon;
        fixedEpsilon.accumulationRate = accumulation_epsilon;
        xGShearDominated = groundingLineSweep( 'shearDominated', ...
            'epsilon',epsilonValues, 'Lambda',LambdaColour, fixedEpsilon);
        figureHandle = figure( 'Color','w', 'Units','centimeters', ...
            'Position',[2,2, ...
            exportStyle.figureWidthCm,exportStyle.figureHeightCm]);
        ax = axes(figureHandle);
        loglog(ax,epsilonValues,xGShearDominated);
        formatSweepAxes(ax,'$\epsilon$',exportStyle);
        colourGradient = applyColourGradient( ...
            ax,LambdaColour,'$\Lambda$',exportStyle);
        epsilonExponent = expectedPowerLawExponent('epsilon',nCurrent);
        xGExtensional = [];
        if extensional.enabled
            xGExtensional = groundingLineSweep( 'extensional', ...
                'epsilon',epsilonValues, 'Lambda',LambdaColour, fixedEpsilon);
            plotExtensionalFamily( ax,epsilonValues,xGExtensional, ...
                colourGradient,extensional.lineWidth,extensional.lineStyle);
            % No separate extensional-model triangle is drawn on the
            % epsilon sweep because both models have the same theoretical
            % power law x_G \propto epsilon^{-1}, so a single triangle is
            % sufficient.
        end
        setSweepAxisLimits( ...
            ax,epsilonValues,xGShearDominated, ...
            xGExtensional,exportStyle);
        shearPlacement = selectTrianglePlacement( ...
            trianglePlacement,'epsilon','shearDominated',nCurrent,nPlotValues);
        addPowerLawTriangle( ...
            ax,epsilonValues,xGShearDominated,epsilonExponent, ...
            shearPlacement, triangleStyle.shearDominated,exportStyle);
        outputFilename = sprintf( 'grounding_line_epsilon_sweep_n%g.png',nCurrent);
        exportgraphics( figureHandle,resultsPath(outputFilename), ...
            'Resolution',exportResolution);
    end

    %% LOCAL SWEEP / PLOTTING FUNCTIONS
    %  The helper functions below select parameter values, call the
    %  grounding-line solvers, and format the resulting parameter sweeps.

    function xGValues = groundingLineSweep( modelName, ...
        xParameterName,xValues, colourParameterName,colourValues, fixed)
        xValues = xValues(:);
        colourValues = colourValues(:).';
        xGValues = NaN(numel(xValues),numel(colourValues));
        for i = 1:numel(xValues)
            for j = 1:numel(colourValues)
                parameters = fixed;
                parameters = setSweepParameter( ...
                    parameters,xParameterName,xValues(i));
                parameters = setSweepParameter( ...
                    parameters,colourParameterName,colourValues(j));
                try
                    switch modelName
                        case 'shearDominated'
                            if parameters.n == 1
                                [xG,~,~,~] = ...
                    findgroundinglineSheardominatednewtonian( ...
                                    parameters.Lambda, parameters.beta, ...
                                    parameters.accumulationRate, ...
                                    parameters.epsilon);
                            else
                                [xG,~,~,~] = ...
                    findgroundinglineSheardominated( ...
                                    parameters.Lambda, parameters.beta, ...
                                    parameters.accumulationRate, ...
                                    parameters.n, parameters.epsilon);
                            end
                        case 'extensional'
                            [xG,~,~,~] = findgroundinglineExtensional( ...
                                parameters.Lambda, parameters.beta, ...
                                parameters.accumulationRate, parameters.n, ...
                                parameters.epsilon);
                        otherwise
                            error( ...
                                'Unknown grounding-line model key: %s.', ...
                                modelName);
                    end
                    xGValues(i,j) = xG;
                catch ME
                    warning( ...
                        ['Grounding-line solve failed for %s model, ', ...
                         '%s = %.4g and %s = %.4g: %s'], modelName, ...
                        xParameterName,xValues(i), ...
                        colourParameterName,colourValues(j), ME.message);
                end
            end
        end
    end

    function parameters = setSweepParameter( ...
        parameters,parameterName,parameterValue)
        switch parameterName
            case 'beta'
                parameters.beta = parameterValue;
            case 'epsilon'
                parameters.epsilon = parameterValue;
            case 'accumulationRate'
                parameters.accumulationRate = parameterValue;
            case 'Lambda'
                parameters.Lambda = parameterValue;
            case 'n'
                parameters.n = parameterValue;
            otherwise
                error('Unknown sweep parameter: %s.',parameterName);
        end
    end

    function plotExtensionalFamily( ...
        ax,xValues,xGValues,colourGradient,lineWidth,lineStyle)
        if size(xGValues,2) ~= size(colourGradient,1)
            error( ...
                ['Number of extensional-model curves does not ', ...
                 'match number of colours.']);
        end
        hold(ax,'on');        for j = 1:size(xGValues,2)
            loglog( ax,xValues,xGValues(:,j), lineStyle, ...
                'Color',colourGradient(j,:), 'LineWidth',lineWidth, ...
                'HandleVisibility','off');
        end
    end

    function formatSweepAxes(ax,xLabelText,exportStyle)
        xlabel(ax,xLabelText, ...
            'Interpreter','latex', ...
            'FontSize',exportStyle.labelFontSizePt);
        ylabel(ax,'$x_G$', ...
            'Interpreter','latex', ...
            'FontSize',exportStyle.labelFontSizePt);
        grid(ax,'off');
        box(ax,'on');
        ax.TickDir = 'out';
        ax.LineWidth = exportStyle.lineWidth;
        ax.FontSize = exportStyle.axesFontSizePt;
    end

    function colourGradient = applyColourGradient( ...
        ax,varyingValues,colourBarLabel,exportStyle)
        lineHandles = findobj(ax,'Type','Line');
        lineHandles = flipud(lineHandles);
        lightColour = [0.70,0.85,1.00];
        darkColour  = [0.02,0.10,0.35];
        nCurves = numel(varyingValues);
        if numel(lineHandles) ~= nCurves
            error( 'Expected %d plotted curves, but found %d.', ...
                nCurves,numel(lineHandles));
        end
        colourGradient = [ ...
            linspace(lightColour(1),darkColour(1),nCurves).', ...
            linspace(lightColour(2),darkColour(2),nCurves).', ...
            linspace(lightColour(3),darkColour(3),nCurves).' ];
        for j = 1:nCurves
            lineHandles(j).Color = colourGradient(j,:);
            lineHandles(j).LineWidth = 1.20;
        end
        colormap(ax,colourGradient);
        clim(ax,log10([min(varyingValues),max(varyingValues)]));
        cb = colorbar(ax,'eastoutside');
        cb.FontSize = exportStyle.colourBarFontSizePt;
        ax.Units = 'normalized';
        cb.Units = 'normalized';
        ax.PositionConstraint = 'innerposition';
        cb.Position = exportStyle.colourBarPosition;
        ax.Position = exportStyle.axesPosition;
        colourBarTicks = logspace( log10(min(varyingValues)), ...
            log10(max(varyingValues)), 5);
        cb.Ticks = log10(colourBarTicks);
        cb.TickLabels = compose('%.2g',colourBarTicks);
        cb.TickDirection = 'out';
        cb.Label.String = '';
        title(cb,colourBarLabel, 'Interpreter','latex', ...
            'FontSize',exportStyle.labelFontSizePt, ...
            'HorizontalAlignment','center');
    end

    function exponent = expectedPowerLawExponent(parameterName,nRheology)
        switch parameterName
            case 'Lambda'
                exponent = -1/(nRheology + 1);
            case 'beta'
                exponent = -(2*nRheology + 1)/(nRheology + 1);
            case 'accumulationRate'
                exponent = 1/(nRheology + 1);
            case 'epsilon'
                exponent = -1;
            otherwise
                error( 'No expected power law has been specified for %s.', ...
                    parameterName);
        end
    end

    function setSweepAxisLimits( ...
        ax,xValues,primaryYValues, ...
        secondaryYValues,exportStyle)
        xValues = xValues(:);
        xValues = xValues(isfinite(xValues) & xValues > 0);
        if ~isempty(xValues)
            xlim(ax,[min(xValues), max(xValues)]);
        end
        allYValues = [primaryYValues(:); secondaryYValues(:)];
        allYValues = allYValues(isfinite(allYValues) & allYValues > 0);
        if ~isempty(allYValues)
            yMinimum = min(allYValues);
            yMaximum = max(allYValues);
            logYMinimum = log10(yMinimum);
            logYMaximum = log10(yMaximum);
            logYRange = logYMaximum - logYMinimum;
            if logYRange <= 0
                logYRange = 1;
            end
            paddedLogYMinimum = ...
                logYMinimum ...
                - exportStyle.lowerYLimitPadding*logYRange;
            if isfield(exportStyle,'upperYLimitPadding')
                paddedLogYMaximum = ...
                    logYMaximum ...
                    + exportStyle.upperYLimitPadding*logYRange;
                yUpperLimit = 10^paddedLogYMaximum;
            else
                yUpperLimit = yMaximum;
            end
            ylim(ax,[10^paddedLogYMinimum, yUpperLimit]);
        end
    end

    function placement = selectTrianglePlacement( ...
        trianglePlacement,parameterName,modelName,nCurrent,nPlotValues)
        index = find(nPlotValues == nCurrent,1,'first');
        if isempty(index)
            error('Could not find n = %g in nPlotValues.',nCurrent);
        end
        rawPlacement = trianglePlacement.(parameterName).(modelName);
        placement.xFraction = rawPlacement.xFraction(index,:);
        if iscell(rawPlacement.relativeTo)
            placement.relativeTo = rawPlacement.relativeTo{index};
        else
            placement.relativeTo = rawPlacement.relativeTo;
        end
        placement.gapFraction = rawPlacement.gapFraction(index);
        if isfield(rawPlacement,'verticalAnchorFraction')
            placement.verticalAnchorFraction = ...
                rawPlacement.verticalAnchorFraction(index);
        end
        if isfield(rawPlacement,'sizeScale')
            placement.sizeScale = rawPlacement.sizeScale(index);
        end
        if iscell(rawPlacement.orientation)
            placement.orientation = rawPlacement.orientation{index};
        else
            placement.orientation = rawPlacement.orientation;
        end
    end

    function exponent = configuredExtensionalExponent( ...
        extensionalExponents,parameterName,nCurrent,nPlotValues)
        index = find(nPlotValues == nCurrent,1,'first');
        if isempty(index)
            exponent = NaN;
            return;
        end
        exponent = extensionalExponents.(parameterName)(index);
    end

    function addPowerLawTriangle( ...
        ax,xValues,yValues,exponent,placement,triangleStyle,exportStyle)
    %ADDPOWERLAWTRIANGLE
    % Draw a theoretical log-log slope triangle.
    % The labels are placed outside the triangle while remaining within
    % the plot limits.
        drawnow;
        xLimits = xlim(ax);
        yLimits = ylim(ax);
        logXLimits = log10(xLimits);
        logYLimits = log10(yLimits);
        logXRange = diff(logXLimits);
        logYRange = diff(logYLimits);

    %% Horizontal triangle position
        logX1 = logXLimits(1) + placement.xFraction(1)*logXRange;
        logX2 = logXLimits(1) + placement.xFraction(2)*logXRange;
        xMarginLeft  = 0.14*logXRange;
        xMarginRight = 0.08*logXRange;
        if logX1 < logXLimits(1) + xMarginLeft
            shift = logXLimits(1) + xMarginLeft - logX1;
            logX1 = logX1 + shift;
            logX2 = logX2 + shift;
        end
        if logX2 > logXLimits(2) - xMarginRight
            shift = logX2 - (logXLimits(2)-xMarginRight);
            logX1 = logX1 - shift;
            logX2 = logX2 - shift;
        end

    %% Triangle rise/run
        [riseInteger,runInteger] = theoreticalIntegerRatio(abs(exponent));
        if ~isfield(placement,'sizeScale')
            placement.sizeScale = 1.0;
        end
        triangleLogXMid = 0.5*(logX1 + logX2);
        triangleWidth = placement.sizeScale*(logX2-logX1);
        logX1 = triangleLogXMid - 0.5*triangleWidth;
        logX2 = triangleLogXMid + 0.5*triangleWidth;
        if logX1 < logXLimits(1) + xMarginLeft
            horizontalCorrection = (logXLimits(1) + xMarginLeft) - logX1;
            logX1 = logX1 + horizontalCorrection;
            logX2 = logX2 + horizontalCorrection;
        end
        if logX2 > logXLimits(2) - xMarginRight
            horizontalCorrection = logX2 - (logXLimits(2) - xMarginRight);
            logX1 = logX1 - horizontalCorrection;
            logX2 = logX2 - horizontalCorrection;
        end
        triangleWidth = logX2 - logX1;
        triangleHeight = abs(exponent)*triangleWidth;

    %% Local plotted curves
        x = xValues(:);
        localRows = x >= 10^logX1 & x <= 10^logX2;
        localY = yValues(localRows,:);
        localY = localY( isfinite(localY) & localY > 0);
        if isempty(localY)
            warning( ['No valid curve values found for ', ...
                 'theoretical triangle placement.']);
            return;
        end

    %% Vertical triangle position
        localLogY = log10(localY);
        localLogYMinimum = min(localLogY);
        localLogYMaximum = max(localLogY);
        if isfield(placement,'verticalAnchorFraction')
            % Prescribed placement used for the accumulation panels.
            logYMid = logYLimits(1) ...
                + placement.verticalAnchorFraction*logYRange;
            if exponent >= 0
                logY1 = logYMid - 0.5*triangleHeight;
                logY2 = logYMid + 0.5*triangleHeight;
            else
                logY1 = logYMid + 0.5*triangleHeight;
                logY2 = logYMid - 0.5*triangleHeight;
            end
        else
            if ~isfield(placement,'gapFraction')
                placement.gapFraction = 0.05;
            end
            verticalGap = placement.gapFraction*logYRange;
            if ~isfield(placement,'relativeTo')
                placement.relativeTo = 'below';
            end
            switch placement.relativeTo
                case 'below'
                    triangleLogTop = localLogYMinimum - verticalGap;
                    if exponent >= 0
                        logY2 = triangleLogTop;
                        logY1 = logY2 - triangleHeight;
                    else
                        logY1 = triangleLogTop;
                        logY2 = logY1 - triangleHeight;
                    end
                case 'above'
                    triangleLogBottom = localLogYMaximum + verticalGap;
                    if exponent >= 0
                        logY1 = triangleLogBottom;
                        logY2 = logY1 + triangleHeight;
                    else
                        logY2 = triangleLogBottom;
                        logY1 = logY2 + triangleHeight;
                    end
                case 'below-upper'
                    triangleLogTop = localLogYMaximum - verticalGap;
                    if exponent >= 0
                        logY2 = triangleLogTop;
                        logY1 = logY2 - triangleHeight;
                    else
                        logY1 = triangleLogTop;
                        logY2 = logY1 - triangleHeight;
                    end
                otherwise
                    error( ...
                        'Unknown relative triangle placement: %s.', ...
                        placement.relativeTo);
            end
        end
        currentLogYLimits = log10(ylim(ax));
        labelSafeMargin = 0.10*diff(currentLogYLimits);
        triangleLogYMinimum = min(logY1,logY2);
        triangleLogYMaximum = max(logY1,logY2);
        if triangleLogYMinimum < currentLogYLimits(1) + labelSafeMargin
            verticalCorrection = (currentLogYLimits(1) + labelSafeMargin) ...
                - triangleLogYMinimum;
            logY1 = logY1 + verticalCorrection;
            logY2 = logY2 + verticalCorrection;
        end
        if triangleLogYMaximum > currentLogYLimits(2) - labelSafeMargin
            verticalCorrection = triangleLogYMaximum ...
                - (currentLogYLimits(2) - labelSafeMargin);
            logY1 = logY1 - verticalCorrection;
            logY2 = logY2 - verticalCorrection;
        end

    %% Draw triangle
        x1 = 10^logX1;
        x2 = 10^logX2;
        y1 = 10^logY1;
        y2 = 10^logY2;
        hold(ax,'on');
        plot(ax,[x1,x2],[y1,y2], 'Color',triangleStyle.Color, ...
            'LineStyle',triangleStyle.LineStyle, ...
            'LineWidth',triangleStyle.LineWidth); % hypotenuse
        if ~isfield(placement,'orientation')
            placement.orientation = 'standard';
        end
        switch placement.orientation
            case 'standard'
                % Right angle at (x1,y2): this is the orientation used in
                % the original script.
                verticalSideX = x1;
                horizontalSideY = y2;
            case 'opposite'
                % Mirror the triangle across its hypotenuse by moving the
                % right angle to (x2,y1).
                verticalSideX = x2;
                horizontalSideY = y1;
            otherwise
                error('Unknown triangle orientation: %s.',placement.orientation);
        end
        plot(ax,[verticalSideX,verticalSideX],[y1,y2], ...
            'Color',triangleStyle.Color, ...
            'LineStyle',triangleStyle.LineStyle, ...
            'LineWidth',triangleStyle.LineWidth); % vertical side
        plot(ax,[x1,x2],[horizontalSideY,horizontalSideY], ...
            'Color',triangleStyle.Color, ...
            'LineStyle',triangleStyle.LineStyle, ...
            'LineWidth',triangleStyle.LineWidth); % horizontal side

    %% Labels outside the triangle, but within plot limits
        currentLogXLimits = log10(xlim(ax));
        currentLogYLimits = log10(ylim(ax));
        currentLogYRange  = diff(currentLogYLimits);
        xLabelMargin = 0.03*logXRange;
        yLabelMargin = 0.045*currentLogYRange;
        logXMid = 0.5*(logX1+logX2);
        logYMid = 0.5*(logY1+logY2);
        if strcmp(placement.orientation,'standard')
            verticalLabelX = logX1 - 0.018*logXRange;
            verticalHA = 'right';
            horizontalSideLogY = logY2;
            if exponent > 0
                horizontalLabelY = horizontalSideLogY + 0.015*currentLogYRange;
                horizontalVA = 'bottom';
            else
                horizontalLabelY = horizontalSideLogY - 0.015*currentLogYRange;
                horizontalVA = 'top';
            end
        else
            verticalLabelX = logX2 + 0.018*logXRange;
            verticalHA = 'left';
            horizontalSideLogY = logY1;
            if exponent > 0
                horizontalLabelY = horizontalSideLogY - 0.015*currentLogYRange;
                horizontalVA = 'top';
            else
                horizontalLabelY = horizontalSideLogY + 0.015*currentLogYRange;
                horizontalVA = 'bottom';
            end
        end
        verticalLabelX = min( ...
            max(verticalLabelX,currentLogXLimits(1) + xLabelMargin), ...
            currentLogXLimits(2) - xLabelMargin);
        verticalLabelY = min( ...
            max(logYMid, currentLogYLimits(1) + yLabelMargin), ...
            currentLogYLimits(2) - yLabelMargin);
        horizontalLabelY = min( ...
            max(horizontalLabelY, currentLogYLimits(1) + yLabelMargin), ...
            currentLogYLimits(2) - yLabelMargin);
        horizontalLabelX = min( ...
            max(logXMid, currentLogXLimits(1) + xLabelMargin), ...
            currentLogXLimits(2) - xLabelMargin);

    %% Draw labels
        text( ax, 10^horizontalLabelX, 10^horizontalLabelY, ...
            sprintf('%d',runInteger), 'HorizontalAlignment','center', ...
            'VerticalAlignment',horizontalVA, ...
            'FontSize',exportStyle.triangleFontSizePt, ...
            'Color',triangleStyle.Color, 'Clipping','on');
        text( ax, 10^verticalLabelX, 10^verticalLabelY, ...
            sprintf('%d',riseInteger), 'HorizontalAlignment',verticalHA, ...
            'VerticalAlignment','middle', ...
            'FontSize',exportStyle.triangleFontSizePt, ...
            'Color',triangleStyle.Color, 'Clipping','on');
    end

    function [numerator,denominator] = theoreticalIntegerRatio(value)
        if ~isfinite(value) || value <= 0
            numerator = 0;
            denominator = 1;
            return;
        end
        [numerator,denominator] = rat(value,1e-10);
        commonFactor = gcd(abs(numerator),abs(denominator));
        numerator = abs(numerator/commonFactor);
        denominator = abs(denominator/commonFactor);
    end
end

%% CORE MODEL FUNCTIONS
% The functions below are shared by all four driver routines.

function [xG,H,h,b] = ...
findgroundinglineSheardominatednewtonian( ...
Lambda,beta,Accumulation,epsilon)
%%

x = @(H) H/(beta + epsilon);

hxFun = @(H) -( ...
    (1 + Accumulation*x(H)) ./ ...
    (H.^3 + Lambda.*H.^2));

F = @(H) ...
    (3/2).*H.^2 ...
    - hxFun(H).^2 ...
    .* ((3/2).*H.^2 + Lambda.*H) ...
    - beta.*Lambda.*hxFun(H).*H ...
    - Accumulation;

%% ------------------------------------------
% Find a positive root only
% -------------------------------------------

Hmin = 1e-10;
Hmax = 1e8;
Hguess = 1;

logHGrid = linspace( ...
    log(Hmin),log(Hmax),2000);

FGrid = arrayfun( ...
    @(logH) F(exp(logH)),logHGrid);

valid = isfinite(FGrid);

idx = find( ...
    valid(1:end-1) ...
    & valid(2:end) ...
    & FGrid(1:end-1).*FGrid(2:end) <= 0);

if isempty(idx)
    error( ...
        ['No positive grounding-line root found ', ...
         'for the Newtonian shear-dominated ', ...
         'ice-sheet model.']);
end

% If more than one positive root exists,
% choose the one whose bracket is closest
% to H = Hguess.
logHmid = 0.5.*( ...
    logHGrid(idx) + logHGrid(idx+1));

[~,j] = min(abs(logHmid - log(Hguess)));
k = idx(j);

logH = fzero( ...
    @(logH) F(exp(logH)), ...
    [logHGrid(k),logHGrid(k+1)]);

H = exp(logH);

%% Grounding-line geometry

xG = H/(beta + epsilon);
h = epsilon*xG;
b = -beta*xG;
end

function [xG,H,h,b] = ...
findgroundinglineSheardominated( ...
Lambda,beta,Accumulation,n,epsilon)
%%

x = @(H) H/(beta + epsilon);

slipCoeff = Lambda;

hxFun = @(H) -( ...
    (1 + Accumulation*x(H)) ./ ...
    (H.^(n+2) ...
    + slipCoeff.*H.^(n+1)) ...
    ).^(1/n);

F = @(H) ...
    Accumulation ...
    + hxFun(H).*abs(hxFun(H)).^(n-1) ...
    .* ( ...
        (n+2)/(n+1).*H.^(n+1).*hxFun(H) ...
        + slipCoeff.*H.^n ...
        .* (hxFun(H) + beta)) ...
    - (n+2)/(n+1).*H.^(n+1);

%% ------------------------------------------
% Find a positive root only
% -------------------------------------------

Hmin = 1e-10;
Hmax = 1e8;
Hguess = 1;

logHGrid = linspace( ...
    log(Hmin),log(Hmax),2000);

FGrid = arrayfun( ...
    @(logH) F(exp(logH)),logHGrid);

valid = isfinite(FGrid);

idx = find( ...
    valid(1:end-1) ...
    & valid(2:end) ...
    & FGrid(1:end-1).*FGrid(2:end) <= 0);

if isempty(idx)
    error( ...
        ['No positive grounding-line root found ', ...
         'for the shear-dominated ice-sheet model.']);
end

logHmid = 0.5.*( ...
    logHGrid(idx) + logHGrid(idx+1));

[~,j] = min(abs(logHmid - log(Hguess)));
k = idx(j);

logH = fzero( ...
    @(logH) F(exp(logH)), ...
    [logHGrid(k),logHGrid(k+1)]);

H = exp(logH);

% Apply grounding-line geometry from
% flotation and b = -beta*x.
xG = H/(beta + epsilon);
h = epsilon*xG;
b = -beta*xG;
end

function [xG,H,h,b] = ...
findgroundinglineExtensional( ...
Lambda,beta,Accumulation,n,epsilon)

if Lambda <= 0
    error( ...
        ['Lambda must be positive for ', ...
         'extensional ice-sheet model.']);
end

p = (n^2 + 3*n + 1)/(n + 1);

K = ((n+2)/(n+1))^(n/(n+1)) ...
    * Lambda^(1/(n+1));

F = @(H) ...
    1 ...
    + Accumulation*H/(beta + epsilon) ...
    - K*H.^p;

%% ------------------------------------------
% Find a positive root only
% -------------------------------------------

Hmin = 1e-10;
Hmax = 1e8;

logHGrid = linspace( ...
    log(Hmin),log(Hmax),2000);

FGrid = arrayfun( ...
    @(logH) F(exp(logH)),logHGrid);

valid = isfinite(FGrid);

idx = find( ...
    valid(1:end-1) ...
    & valid(2:end) ...
    & FGrid(1:end-1).*FGrid(2:end) <= 0, ...
    1,'first');

if isempty(idx)
    error( ...
        ['No positive grounding-line root found ', ...
         'for the extensional ice-sheet model.']);
end

logH = fzero( ...
    @(logH) F(exp(logH)), ...
    [logHGrid(idx),logHGrid(idx+1)]);

H = exp(logH);

%% Grounding-line geometry

xG = H/(beta + epsilon);
h = epsilon*xG;
b = -beta*xG;
end

function [ ...
xsheet,xshelf,hsheet,hshelf, ...
bsheet,bshelf,hxsheet,xG,HG,hG,bG] = ...
Sheardominatedicesheetnewtonian( ...
Lambda,beta,Accumulation, ...
epsilon,shelfextent)
%%
% Sheardominatedicesheetnewtonian solves the
% governing equations for the Newtonian
% shear-dominated ice-sheet model.
%
% It returns h and b for both the sheet and shelf,
% the x values in both regions, and the grounding-
% line quantities xG, HG, hG and bG.

% Find grounding-line position.
[xG,HG,hG,bG] = ...
findgroundinglineSheardominatednewtonian( ...
Lambda,beta,Accumulation,epsilon);

% Solve sheet equations upstream from x_G to 0.
[xsheet,Hsheet] = ode45( ...
@(x,H) sheetequationSheardominatednewtonian( ...
x,H,Lambda,beta,Accumulation), ...
[xG,0],HG);

% Apply b = -Ax and h = H+b.
bsheet = -beta*xsheet;
hsheet = Hsheet + bsheet;

% h_x = H_x - beta, used to determine velocity.
hxsheet = -( ...
(1 + Accumulation*xsheet) ./ ...
(Hsheet.^3 + Lambda.*Hsheet.^2));

% Solve shelf equations downstream from x_G.
[xshelf,Hshelf] = ode45( ...
@(x,H) shelfequationSheardominatednewtonian( ...
x,H,Accumulation), ...
[xG,xG*shelfextent],HG);

% Use flotation condition to determine h and b.
flotationh = epsilon/(beta + epsilon);
flotationb = beta/(beta + epsilon);

hshelf = flotationh*Hshelf;
bshelf = -flotationb*Hshelf;

end

function dHdx = ...
sheetequationSheardominatednewtonian( ...
x,H,Lambda,beta,Accumulation)
% Solve the Newtonian shear-dominated
% ice-sheet model sheet-depth equation.

dHdx = ...
    beta ...
    - (1 + Accumulation*x) ./ ...
    (H.^3 + Lambda.*H.^2);
end

function dHdx = ...
shelfequationSheardominatednewtonian( ...
x,H,Accumulation)
% Solve the Newtonian shear-dominated
% ice-sheet model shelf-depth equation.

dHdx = ...
    (Accumulation*H - (3/2).*H.^3) ...
    ./ (1 + Accumulation*x);
end

function [ ...
xsheet,xshelf,hsheet,hshelf, ...
bsheet,bshelf,hxsheet,xG,HG,hG,bG] = ...
Sheardominatedicesheet( ...
Lambda,beta,Accumulation,n, ...
epsilon,shelfextent)
%%
% Sheardominatedicesheet solves the governing
% equations for the non-Newtonian shear-dominated
% ice-sheet model.
%
% It returns h and b for both the sheet and shelf,
% the x values in both regions, and the grounding-
% line quantities xG, HG, hG and bG.

% Find grounding-line position.
[xG,HG,hG,bG] = ...
findgroundinglineSheardominated( ...
Lambda,beta,Accumulation,n,epsilon);

% Solve sheet equations upstream from x_G to 0.
[xsheet,Hsheet] = ode45( ...
@(x,H) sheetequationSheardominated( ...
x,H,Lambda,beta,Accumulation,n), ...
[xG,0],HG);

% Apply b = -Ax and h = H+b.
bsheet = -beta*xsheet;
hsheet = Hsheet + bsheet;

% h_x = H_x - beta, used to determine velocity.
hxsheet = -( ...
(1 + Accumulation*xsheet) ./ ...
(Hsheet.^(n+2) ...
+ Lambda.*Hsheet.^(n+1)) ...
).^(1/n);

% Solve shelf equations downstream from x_G.
[xshelf,Hshelf] = ode45( ...
@(x,H) shelfequationSheardominated( ...
x,H,Accumulation,n), ...
[xG,xG*shelfextent],HG);

% Use flotation condition to determine h and b.
flotationh = epsilon/(beta + epsilon);
flotationb = beta/(beta + epsilon);

hshelf = flotationh*Hshelf;
bshelf = -flotationb*Hshelf;

end

function dHdx = ...
sheetequationSheardominated( ...
x,H,Lambda,beta,Accumulation,n)
% Solve the non-Newtonian shear-dominated
% ice-sheet model sheet-depth equation.

dHdx = ...
    beta ...
    - ( ...
        (1 + Accumulation*x) ./ ...
        (H.^(n+2) ...
        + Lambda.*H.^(n+1)) ...
      ).^(1/n);
end

function dHdx = ...
shelfequationSheardominated( ...
x,H,Accumulation,n)
% Solve the non-Newtonian shear-dominated
% ice-sheet model shelf-depth equation.

dHdx = ...
    (Accumulation*H ...
    - (n+2)/(n+1).*H.^(n+2)) ...
    ./ (1 + Accumulation*x);
end

function [ ...
xsheet,xshelf,hsheet,hshelf, ...
bsheet,bshelf,xG,HG,hG,bG] = ...
Extensionalicesheet( ...
Lambda,beta,Accumulation,n, ...
epsilon,shelfextent)
%%
% Extensionalicesheet solves the governing
% equations for the non-Newtonian extensional
% ice-sheet model.
%
% It returns h and b for both the sheet and shelf,
% the x values in both regions, and the grounding-
% line quantities xG, HG, hG and bG.

% Find grounding-line position.
[xG,HG,hG,bG] = ...
findgroundinglineExtensional( ...
Lambda,beta,Accumulation,n,epsilon);

% Solve sheet equations upstream from x_G to 0.
[xsheet,Hsheet] = ode45( ...
@(x,H) sheetequationExtensional( ...
x,H,Lambda,beta,Accumulation,n), ...
[xG,0],HG);

% Apply b = -Ax and h = H+b.
bsheet = -beta*xsheet;
hsheet = Hsheet + bsheet;

% Solve shelf equations downstream from x_G.
[xshelf,Hshelf] = ode45( ...
@(x,H) shelfequationExtensional( ...
x,H,Accumulation,n), ...
[xG,xG*shelfextent],HG);

% Use flotation condition to determine h and b.
flotationh = epsilon/(beta + epsilon);
flotationb = beta/(beta + epsilon);

hshelf = flotationh*Hshelf;
bshelf = -flotationb*Hshelf;

end

function dHdx = ...
sheetequationExtensional( ...
x,H,Lambda,beta,Accumulation,n)
% Solve the non-Newtonian extensional
% ice-sheet model sheet-depth equation.

dHdx = ...
    beta ...
    - ( ...
        (1 + Accumulation*x) ./ ...
        (Lambda.*H.^(n+1)) ...
      ).^(1/n);
end

function dHdx = ...
shelfequationExtensional( ...
x,H,Accumulation,n)
% Solve the extensional ice-sheet
% shelf-depth equation.

dHdx = ...
    (Accumulation*H ...
    - (n+2)/(n+1).*H.^(n+2)) ...
    ./ (1 + Accumulation*x);
end

function [ xsheet,xshelf,hsheet,hshelf, bsheet,bshelf,hxsheet,xG, ...
    HG,hG,bG,Accumulation] = SheardominatedicesheetConstantShelf( ...
        Lambda,beta,n, epsilon,shelfextent)
%%
% SheardominatedicesheetConstantShelf solves the governing equations
% for the shear-dominated ice-sheet model subject to constant shelf depth.
% It returns h and b for both the sheet and shelf,
% the x values in both regions, the sheet surface gradient,
% the grounding-line quantities xG, HG, hG and bG,
% and the accumulation rate required for constant shelf depth.

% Find grounding-line position and accumulation rate.
[xG,HG,hG,bG,Accumulation] = ...
    findgroundinglineConstantShelfSheardominated( Lambda,beta,n,epsilon);

% Solve sheet equations upstream from x_G to 0.
[xsheet,Hsheet] = ode45( @(x,H) sheetequationSheardominated( ...
        x,H,Lambda,beta,Accumulation,n), [xG,0], HG);

% Apply b = -beta*x and h = H+b.
bsheet = -beta*xsheet;
hsheet = Hsheet + bsheet;

% Surface gradient h_x.
hxsheet = -( (1 + Accumulation*xsheet) ./ ( Hsheet.^(n+2) ...
        + Lambda.*Hsheet.^(n+1) ) ).^(1/n);

% Solve shelf equations downstream from x_G.
[xshelf,Hshelf] = ode45( @(x,H) shelfequationSheardominated( ...
        x,H,Accumulation,n), [xG,xG*shelfextent], HG);

% Use flotation condition to determine shelf h and b.
flotationh = epsilon/(beta + epsilon);
flotationb = beta/(beta + epsilon);

hshelf = flotationh*Hshelf;
bshelf = -flotationb*Hshelf;

end

function [xG,H,h,b,Accumulation] = ...
    findgroundinglineConstantShelfSheardominated( Lambda,beta,n,epsilon)
%%
% Determine the grounding-line depth for a constant-depth shelf.

    x = @(H) H/(beta + epsilon);

    % Constant-shelf-depth accumulation condition.
    % Accumulation = (n+2)/(n+1) H_G^(n+1).
    accumulationFun = @(H) (n+2)/(n+1) .* H.^(n+1);

    % Surface gradient at the grounding line.
    hxFun = @(H) -( ( 1 + accumulationFun(H).*x(H) ) ./ ( H.^(n+2) ...
            + Lambda.*H.^(n+1) ) ).^(1/n);

    % Full grounding-line stress condition.
    F = @(H) accumulationFun(H) + hxFun(H) .* abs(hxFun(H)).^(n-1) ...
        .* ( (n+2)/(n+1) .* H.^(n+1) .* hxFun(H) + Lambda .* H.^n ...
            .* (hxFun(H) + beta) ) - (n+2)/(n+1) .* H.^(n+1);

    %% Find a positive grounding-line thickness

    Hmin = 1e-10;
    Hmax = 1e8;
    Hguess = 1;

    logHGrid = linspace( log(Hmin), log(Hmax), 2000);

    FGrid = arrayfun( @(logH) F(exp(logH)), logHGrid);

    valid = isfinite(FGrid);

    idx = find( valid(1:end-1) & valid(2:end) & FGrid(1:end-1) ...
        .* FGrid(2:end) <= 0);

    if isempty(idx)
        error( ['No positive constant-shelf shear-dominated ', ...
             'ice-sheet grounding-line root found.']);
    end

    logHmid = 0.5.*( logHGrid(idx) + logHGrid(idx+1));

    [~,j] = min( abs(logHmid - log(Hguess)));

    k = idx(j);

    logH = fzero( @(logH) F(exp(logH)), [ logHGrid(k), logHGrid(k+1) ]);

    H = exp(logH);

    %% Grounding-line geometry and accumulation

    Accumulation = accumulationFun(H);

    xG = H/(beta + epsilon);
    h = epsilon*xG;
    b = -beta*xG;

end

function [ xsheet,xshelf,hsheet,hshelf, ...
    bsheet,bshelf,xG,HG,hG,bG,Accumulation] = ...
    ExtensionalicesheetConstantShelf( Lambda,beta,n, epsilon,shelfextent)
%%
% ExtensionalicesheetConstantShelf solves the governing equations
% for the extensional ice-sheet model subject to constant shelf depth.

% Determine the grounding-line position and accumulation rate.
[xG,HG,hG,bG,Accumulation] = ...
    findgroundinglineConstantShelfExtensional( Lambda,beta,n,epsilon);

% Solve sheet equations upstream from x_G to 0.
[xsheet,Hsheet] = ode45( @(x,H) sheetequationExtensional( ...
        x,H,Lambda,beta,Accumulation,n), [xG,0], HG);

% Apply b = -beta*x and h = H+b.
bsheet = -beta*xsheet;
hsheet = Hsheet + bsheet;

% Solve shelf equations downstream from x_G.
[xshelf,Hshelf] = ode45( @(x,H) shelfequationExtensional( ...
        x,H,Accumulation,n), [xG,xG*shelfextent], HG);

% Apply the flotation condition in the shelf.
flotationh = epsilon/(beta + epsilon);
flotationb = beta/(beta + epsilon);

hshelf = flotationh*Hshelf;
bshelf = -flotationb*Hshelf;

end

function [xG,H,h,b,Accumulation] = ...
    findgroundinglineConstantShelfExtensional( Lambda,beta,n,epsilon)

    if Lambda <= 0
        error( ['Lambda must be positive ', ...
             'for the extensional ice-sheet model.']);
    end

    x = @(H) H/(beta + epsilon);

    accumulationFun = @(H) (n+2)/(n+1) .* H.^(n+1);

    p = (n^2 + 3*n + 1) /(n + 1);

    K = ((n+2)/(n+1)) ^(n/(n+1)) .* Lambda ^(1/(n+1));

    % Extensional grounding-line flux condition.
    F = @(H) 1 + accumulationFun(H).*x(H) - K.*H.^p;

    %% Find a positive grounding-line thickness

    Hmin = 1e-10;
    Hmax = 1e8;
    Hguess = 1;

    logHGrid = linspace( log(Hmin), log(Hmax), 2000);

    FGrid = arrayfun( @(logH) F(exp(logH)), logHGrid);

    valid = isfinite(FGrid);

    idx = find( valid(1:end-1) & valid(2:end) & FGrid(1:end-1) ...
        .* FGrid(2:end) <= 0);

    if isempty(idx)
        error( ['No positive constant-shelf extensional ', ...
             'ice-sheet grounding-line root found.']);
    end

    logHmid = 0.5.*( logHGrid(idx) + logHGrid(idx+1));

    [~,j] = min( abs(logHmid - log(Hguess)));

    k = idx(j);

    logH = fzero( @(logH) F(exp(logH)), [ logHGrid(k), logHGrid(k+1) ]);

    H = exp(logH);

    %% Grounding-line geometry and accumulation

    Accumulation = accumulationFun(H);

    xG = H/(beta + epsilon);    h = epsilon*xG;
    b = -beta*xG;

end

function u = ...
velocityprofileSheardominatedsheetnewtonian( ...
eta,H,hx,Lambda)
%VELOCITYPROFILESHEARDOMINATEDSHEETNEWTONIAN
%
% Nondimensional Newtonian shear-dominated
% ice-sheet model velocity profile.
%
% eta = 0 at the bed
% eta = 1 at the surface
%
% hx = dh/dx < 0

s = -hx;
profileCoeff = 3/2;

% Basal sliding velocity.
ub = Lambda.*H.*s;

% Full velocity profile.
u = ...
    ub ...
    + profileCoeff.*H.^2.*s ...
    .* (1 - (1-eta).^2);
end

function u = ...
velocityprofileSheardominatedsheet( ...
eta,H,hx,Lambda,n)
%VELOCITYPROFILESHEARDOMINATEDSHEET
%
% Nondimensional shear-dominated
% ice-sheet model velocity profile.
%
% eta = 0 at the bed
% eta = 1 at the surface
%
% hx = dh/dx < 0

s = -hx;
profileCoeff = (n+2)/(n+1);

% Basal sliding velocity.
ub = Lambda.*H.^n.*s.^n;

% Full velocity profile.
u = ...
    ub ...
    + profileCoeff.*H.^(n+1).*s.^n ...
    .* (1 - (1-eta).^(n+1));
end

function u = ...
velocityprofileExtensionalsheet( ...
eta,x,H,Accumulation)
%VELOCITYPROFILEEXTENSIONALSHEET
%
% Extensional ice-sheet model outer solution.
% Velocity is independent of depth.
%
% eta is included so that this function
% has the same style as the shear-dominated
% ice-sheet model velocity-profile function.

q = 1 + Accumulation*x;
u0 = q./H;

u = u0.*ones(size(eta));
end

function u = ...
velocityprofileshelf( ...
eta,x,H,Accumulation)
%VELOCITYPROFILESHELF
%
% Depth-independent shelf velocity.

q = 1 + Accumulation*x;
u0 = q./H;

u = u0.*ones(size(eta));
end
%% OUTPUT PATH HELPER
function outputPath = resultsPath(filename)
%RESULTSPATH Return an output path in the results folder beside this file.
    driverDirectory = fileparts(mfilename('fullpath'));
    outputPath = fullfile(driverDirectory, 'results', filename);
end
