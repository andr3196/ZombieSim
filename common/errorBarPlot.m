function res = errorBarPlot(figureStruct) %x,y,yfit,xerror, yerror, xlab,ylab, filename, newfig)
%ERRORBARPLOT 2-D figure plotting
%   ERRORBARPLOT(F) plots multiple datasets with vertical and horizontal errorbars, any
%   number of fitted curves, sets axis labels and ... F is a MATLAB
%   STRUCTURE
%   Allowed fields:
%   'x' values of the independent variable (matrix of arbitrary
%   orientation)
%   'y' values of the dependent variable  (matrix of same size as x)
%   'xfit' values of the independent variable (matrix of arbitrary orientation)
%   'theta' values of the independent variable in polar coordinates (matrix)
%   'rho' values of the dependent variable in polar coordinates (matrix with the same orientation as 'theta')
%   'thetaError' errors on the independent variable
%   'rhoError' errors on the dependent variable
%   'thetaFit'
%   'rhoFit'
%   'thetaLimits'
%   'thetaZeroLocation'
%   'polarFitAttributes'
%   'polarAttributes' set 'Color', 'LineWidth', 'MaxHeadSize'... for each
%   data points with or without errorbars !!!!!!!!!!!!!!!!!!!!! NB: NOT
%   IMPLEMENTED YET !!!!!!!!!!!!!!!!!!
%   'yfit' values of fitted curve (matrix of same size as xfit)
%   'xError' errors on the independent variable (matrix of same size as x)
%   'yError' errors on the dependent variable (matrix of the same size as y)
%   'errorBarLineWidth' thickness of the line making up the errorbars (both
%   vertical and horisontal) (integer)
%   'markerSize' in case of no errorbars this value adjusts the size of the
%   symbol being plotted (integer)
%   'fitLineStyle' line colour and marker symbol for the fit. Either single value ('b:') if the same format for all or cell vector of same length as number of fits
%   'fitLineWidth' width of the fitting curves (integer)
%   'xLabel' label for the horizontal axis (string)
%   'yLabel' label for the vertical axis (string)
%   'zLabel' label for the third axis (string)
%   'xlimits'sets the interval on the x-axis ([xmin xmax])
%   'ylimits'sets the interval on the y-axis ([ymin ymax])
%   'legend' legend (cell array)
%   'newFigure' whether to plot data in current or new figure (arbitrary value, maybe boolean)
%   'fontsize' axis fontsize (int)
%   'zoomFrameMethod' whether a zoom of a portion of the data should be
%   shown (either 'Automatic' or 'Manual')
%   'zoomFramex' independent variable in zoom frame (matrix of arbitrary orientation)
%   'zoomFramey' dependent variable in zoom frame (matrix of same size as zoomFramex)
%   'zoomFrameMarkerSize' in case of no errorbars this value adjusts the size of the
%   symbol being plotted in the zoomFrame (integer)
%   'zoomFramexfit' values of the independent variable for the fit in the zoom frame (matrix of arbitrary orientation)
%   'zoomFrameyfit' values of the dependent variable for the fit in the zoom frame (matrix of same size as zoomFramexfit)
%   'zoomFitLineStyle' line colour and marker symbol for the fit in the zoom frame. Either single value ('b:') if the same format for all or cell vector of same length as number of fits
%   'zoomFitLineWidth' width of the fitting curves in zoom frame (integer)
%   'zoomFontSize' font size of zoomframe (integer)
%   'zoomFrameThetaFit'
%   'zoomFrameRhoFit'
%   'zoomFramePolarAttributes'
%   'zoomXLabel' xlabel in zoomFrame (string)
%   'zoomYLabel' ylabel in zoomFrame (string)
%   'zoomZLabel' zlabel in zoomFrame (string)
%   'zoomXlimits' sets the interval on the x-axis in the zoomframe ([xmin xmax])
%   'zoomYlimits' sets the interval on the y-axis in the zoomframe([ymin ymax])
%   'zoomFrameLegend' legend in zoomFrame (cell array)
%   'zoomAxesGrid' whether grid lines should be shown in zoomFrame (boolean)
%   'zoomAxes' axes limits in zoomFrame (vector)
%   'method3D' what type of 3D-plot is wanted (either 'contour' or 'surf')
%   'XY' idependent variables x and y in 3D (matrix of size (n, m, 2, p), where n and m are the dimensions of the X and Y matrices and p is the number of data sets to be plotted)
%   'Z' dependent variable z in 3D (matrix of size (n,m,p))
%   'axesGrid' whether to display axis grid lines (boolean)
%   'axes' axes limits in main frame (vector)
%   'folder' the folder in which the figure should be saved (path as string)
%   'filename' the name of the file in which the figure is saved (string)
%   'filetype' extension of the file
%   'arrows' allows placement of arrows on the figure (double array with
%   the format: [x1 y1 dx1 dy2; x2 y2 ...], where (x,y) are the
%   coordinates of the backend of the arrow and dx and dy are the relative
%   component distances to the head
%   'placeArrows' allows for placement of arrows via a guide (any value)
%   'arrowAttributes' set 'Color', 'LineWidth', 'MaxHeadSize'... for each
%   arrow or all arrows at once (cell array of key-value-pairs either with one row if all arrows should look the same, or number of rows corresponding to the number of rows in 'arrows')
%   'plotAttributes' set 'Color', 'LineWidth', 'MaxHeadSize'... for each
%   data points with or without errorbars (same as arrowAttributes)
%   'fitAttributes' set 'Color', 'LineWidth', 'MaxHeadSize'... for each
%   curve with or without errorbars (same as arrowAttributes)
%   'isSemilogy' set whether to use semilog plot (boolean)
%   'slimOff' (boolean) prevent removal of whitespace 


%-------------- Additional details --------------------------------------
% If data of different dimensions needs to be supported, let all extra
% entries be 'nan'. For instance in order to provide y1 = [1 2 3] and y2 =
% [3 4], let  y = [1 2 3; 3 4 NaN]
%(c) Andreas Springborg 2016

if nargin < 1
    error(['Expected 1 input, but got: ' num2str(nargin)]);
elseif ~isstruct(figureStruct)
    error(['Expected input to be of type ''struct'', but was: ' class(figureStruct)])
end

% Definitions:
global fitLineWidthDefault markerSizeDefault errorBarLineWidthDefault filetypeDefault...
    errorBarMarkerDefault isSemilogDefault;
fitLineWidthDefault = 1;
errorBarLineWidthDefault = 1;
markerSizeDefault = 5;

%%% Errorbars
errorBarMarkerDefault = 'None';

%%% Semilog mode
isSemilogDefault = false;

%%% Files
filetypeDefault = '-deps';

%%% Arrows



% End of definitions

if isfield(figureStruct,'newFigure') % determine whether to start a new figure
    if figureStruct.newFigure
        figure
        if isfield(figureStruct, 'rho')
           polaraxes 
        end
    end
end

if isfield(figureStruct,'figureNumber') % determine whether to start a new figure
    figure(figureStruct.figureNumber)
end

%polarplot(1,1)
hold on

[or, h] = handleErrorBars(figureStruct, 'normal'); % perform actions related to 2D errorbars
[or, h] = handleErrorBars(figureStruct, 'rightaxis'); % perform actions related to 2D errorbars on right axis
h = handle3DPlot(figureStruct, h); % perform actions related to 3D plots
h = handleFittingCurves(figureStruct, or, h, 'normal'); % perform actions related to 2D fitting curves
h = handleFittingCurves(figureStruct, or, h, 'rightaxis'); % perform actions related to 2D fitting curves
[or, h] = handlePolarPlotting(figureStruct, h);
h = handlePolarFittingCurves(figureStruct, or, h, 'normal');
figureStruct.plotHandles = h;


figureStruct = handleAppearance(figureStruct, 'normal'); % performs actions related to setting fontsize, axis labels, and grid lines
figureStruct = handleAppearance(figureStruct, 'rightaxis'); % performs actions related to setting fontsize, axis labels, and grid lines on right axis
figureStruct.axisHandle = gca;
figureStruct.figureHandle = gcf;

handleZoomFrame(figureStruct, or) % perform actions related to placing a zoomed part of the figure in the figure

handleArrows(figureStruct)
handleSlim(figureStruct)
if isfield(figureStruct,'filename')
    if ~isfield(figureStruct, 'filetype')
        figureStruct.filetype = filetypeDefault;
    end
    if isfield(figureStruct, 'folder')
        print(gcf,figureStruct.filetype, [figureStruct.folder filesep figureStruct.filename])
    else
        print(figureStruct.filename,figureStruct.filetype)
    end
end

if nargout > 0
    res = figureStruct;
end



function [or, h] = handleErrorBars(figureStruct, mode)
global markerSizeDefault errorBarLineWidthDefault errorBarMarkerDefault isSemilogDefault
h = [];
or = NaN;
if strcmpi(mode, 'normal')
    yField = 'y';
    xField = 'x';
    yErrorField = 'yError';
    xErrorField = 'xError';
    markerSize = 'markerSize';
    marker = 'marker';
    lineColor = 'errorbarColor';
    errorBarLineWidth = 'errorBarLineWidth';
    semilogMode = 'isSemilogy';
    attributes = 'plotAttributes';
    %yyaxis left
elseif strcmpi(mode, 'rightaxis')
    if isempty(cell2mat(strfind(fieldnames(figureStruct), 'right')))
        return
    end
    yField = 'rightY';
    xField = 'rightX';
    yErrorField = 'rightYError';
    xErrorField = 'rightXError';
    markerSize = 'rightMarkerSize';
    marker = 'rightMarker';
    lineColor = 'rightErrorbarColor';
    errorBarLineWidth = 'rightErrorBarLineWidth';
    semilogMode = 'rightIsSemilogy';
    attributes = 'rightPlotAttributes';
    yyaxis right 
elseif strcmpi(mode, 'zoom')
    yField = 'zoomFramey';
    xField = 'zoomFramex';
    yErrorField = 'zoomFrameyError';
    xErrorField = 'zoomFramexError';
    markerSize = 'zoomFrameMarkerSize';
    marker = 'zoomFrameMarker';
    lineColor = 'zoomFrameErrorBarColor';
    errorBarLineWidth = 'zoomFrameLineWidth';
    semilogMode = 'zoomFrameIsSemilogy';
    attributes = 'zoomFrameAttributes';
else
    return
end
if ~isfield(figureStruct, markerSize)
    figureStruct.(markerSize) = markerSizeDefault;
end
if ~isfield(figureStruct, errorBarLineWidth)
    figureStruct.(errorBarLineWidth) = errorBarLineWidthDefault;
end


if isfield(figureStruct, yField)
    
    y = figureStruct.(yField);
    or = getOrientation(y);
    if or < 0
        y = y.';
    end
    
    if isfield(figureStruct, attributes)
        plotAttributes = figureStruct.(attributes);
        if size(plotAttributes,1)
            plotAttributes = repmat(plotAttributes,size(y,1),1);
        end
    else
        plotAttributes = {};
    end
    
    if ~isfield(figureStruct,semilogMode)
        figureStruct.(semilogMode) = isSemilogDefault;
    end
    
    if ~isfield(figureStruct, marker)
        figureStruct.(marker) = repmat({errorBarMarkerDefault}, 1, size(y,1));
    end
    if ~isfield(figureStruct, lineColor)
        figureStruct.(lineColor) = repmat({rand(1,3)}, 1, size(y,1));
    end
    if isfield(figureStruct, xField)
        x = figureStruct.(xField);
        if or < 0
            x = x.';
        end
        if isfield(figureStruct,yErrorField)
            yError = figureStruct.(yErrorField);
            if or < 0
                yError = yError.';
            end
            for i = 1:size(y,1)
                h(end+1) = errorbar(x(i,:),y(i,:),yError(i,:),'LineStyle','none', 'LineWidth', figureStruct.(errorBarLineWidth), 'Marker', figureStruct.(marker){i}, 'MarkerSize', figureStruct.(markerSize), 'Color', figureStruct.(lineColor){i});
                if ~isempty(plotAttributes)
                    for j = 1:2:size(plotAttributes,2)
                        if ~isempty(plotAttributes{i,j})
                            set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                        end
                    end
                end
            end
        else
            if figureStruct.(semilogMode)
                hold off
                for i = 1: size(y,1)
                    h(end+1) = semilogy(x(i,:),y(i,:),'bx', 'MarkerSize', figureStruct.(markerSize));
                    if i == 1
                        hold on
                    end
                    if ~isempty(plotAttributes)
                        for j = 1:2:size(plotAttributes,2)
                            if ~isempty(plotAttributes{i,j})
                                set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(y,1)
                    h(end+1) = plot(x(i,:),y(i,:),'bx', 'MarkerSize', figureStruct.(markerSize));
                    if ~isempty(plotAttributes)
                        for j = 1:2:size(plotAttributes,2)
                            if ~isempty(plotAttributes{i,j})
                                set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
        end
        if isfield(figureStruct, xErrorField)
            xError = figureStruct.(xErrorField);
            if or < 0
                xError = xError.';
            end
            for i = 1:size(y,1)
                h(end+1) = herrorbar(x(i,:),y(i,:),xError(i,:), 'LineWidth', figureStruct.(errorBarLineWidth), figureStruct.(lineColor){i});
                if ~isempty(plotAttributes)
                    for j = 1:2:size(plotAttributes,2)
                        if ~isempty(plotAttributes{i,j})
                            set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    else
        if isSemilog
            hold off
            for i = 1: size(y,1)
                h(end+1) = semilogy(y(i,y(i,:)),'bx', 'MarkerSize', figureStruct.(markerSize));
                if i == 1
                    hold on
                end
                if ~isempty(plotAttributes)
                    for j = 1:2:size(plotAttributes,2)
                        if ~isempty(plotAttributes{i,j})
                            set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                        end
                    end
                end
            end
        else
            for i = 1:size(y,1)
                h(end+1) = plot(y(i,y(i,:)),'bx', 'MarkerSize', figureStruct.(markerSize));
                if ~isempty(plotAttributes)
                    for j = 1:2:size(plotAttributes,2)
                        if ~isempty(plotAttributes{i,j})
                            set(h(end),plotAttributes{i,j}, plotAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    end
end
if strcmpi(mode, 'rightaxis')
    yyaxis left
end

function h = handleFittingCurves(figureStruct, or, h, mode)
global fitLineWidthDefault
if strcmpi(mode, 'normal')
    xfitField = 'xfit';
    yfitField = 'yfit';
    yfitLineStyle = 'fitLineStyle';
    yfitLineWidth = 'fitLineWidth';
elseif strcmpi(mode, 'rightaxis')
    if isempty(cell2mat(strfind(fieldnames(figureStruct), 'right')))
        return
    end
    xfitField = 'rightXfit';
    yfitField = 'rightYfit';
    yfitLineStyle = 'rightFitLineStyle';
    yfitLineWidth = 'rightFitLineWidth';
    yyaxis right
elseif strcmpi(mode, 'zoom') 
    xfitField = 'zoomFramexfit';
    yfitField = 'zoomFrameyfit';
    yfitLineStyle = 'zoomFitLineStyle';
    yfitLineWidth = 'zoomFitLineWidth';
else
    h = [];
    return
end
if ~ isfield(figureStruct, yfitLineWidth)
    figureStruct.(yfitLineWidth) = fitLineWidthDefault;
end

if  isfield(figureStruct,xfitField)
    xfit = figureStruct.(xfitField);
    if or < 0
        xfit = xfit.';
    end
    if isfield(figureStruct, 'fitAttributes')
        fitAttributes = figureStruct.('fitAttributes');
        if size(fitAttributes,1)
            fitAttributes = repmat(fitAttributes,size(xfit,1),1);
        end
    else
        fitAttributes = {};
    end
    
    if isfield(figureStruct, yfitField)
        yfit = figureStruct.(yfitField);
        if or < 0
            yfit = yfit.';
        end
        if isfield(figureStruct,yfitLineStyle)
            lineStyle = figureStruct.(yfitLineStyle);
            if iscell(lineStyle)
                for i = 1:size(yfit,1)
                    h(end+1) = plot(xfit(i,:),yfit(i,:),lineStyle{i}, 'LineWidth', figureStruct.(yfitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:size(fitAttributes,2)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(yfit,1)
                    h(end+1) = plot(xfit(i,:),yfit(i,:),lineStyle, 'LineWidth', figureStruct.(yfitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:size(fitAttributes,2)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
            
        else
            for i = 1:size(yfit,1)
                h(end+1) = plot(xfit(i,:),yfit(i,:),'r-', 'LineWidth', figureStruct.(yfitLineWidth));
                if ~isempty(fitAttributes)
                    for j = 1:2:size(fitAttributes,2)
                        if ~isempty(fitAttributes{i,j})
                            set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    end
else
    if isfield(figureStruct, yfitField)
        yfit = figureStruct.(yfitField);
        if or < 0
            yfit = yfit.';
        end
        if isfield(figureStruct,yfitLineStyle)
            lineStyle = figureStruct.(yfitLineStyle);
            if iscell(lineStyle)
                for i = 1:size(yfit,1)
                    h(end+1) = plot(yfit(i,:),lineStyle{i}, 'LineWidth', figureStruct.(yfitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(yfit,1)
                    h(end+1) = plot(yfit(i,:),lineStyle, 'LineWidth', figureStruct.(yfitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
            
        else
            for i = 1:size(yfit,1)
                h(end+1) = plot(yfit(i,:),'r-', 'LineWidth', figureStruct.(yfitLineWidth));
                if ~isempty(fitAttributes)
                    for j = 1:2:length(fitAttributes)
                        if ~isempty(fitAttributes{i,j})
                            set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    end
end
if strcmpi(mode, 'rightaxis')
    yyaxis left
end

function [or, h] = handlePolarPlotting(figureStruct, h)
if isfield(figureStruct, 'polarAttributes')
    polarAttributes = figureStruct.(attributes);
    if size(polarAttributes,1)
        polarAttributes = repmat(polarAttributes,size(y,1),1);
    end
else
    polarAttributes = {};
end
or = 0;
if isfield(figureStruct, 'rho')
    rho = figureStruct.rho;
    or = getOrientation(rho);
    if isfield(figureStruct, 'theta')
        theta = figureStruct.theta;
        if isfield(figureStruct, 'rhoError')
            rhoError = figureStruct.rhoError;
            if isfield(figureStruct, 'thetaError')
                thetaError = figureStruct.thetaError;
                for i = 1:size(theta,1)
                    perrorbar(theta(i,:), rho(i,:),thetaError(i,:), rhoError(i,:));
                    if ~isempty(polarAttributes)
                        for j = 1:2:size(polarAttributes,2)
                            if ~isempty(polarAttributes{i,j})
                                set(h(end),polarAttributes{i,j}, polarAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(theta,1)
                    perrorbar(theta(i,:), rho(i,:),[], rhoError(i,:));
                    if ~isempty(polarAttributes)
                        for j = 1:2:size(polarAttributes,2)
                            if ~isempty(polarAttributes{i,j})
                                set(h(end),polarAttributes{i,j}, polarAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
        else
            for i = 1:size(theta,1)
                h(end+1) =  polarplot(theta(i,:), rho(i,:));
                if ~isempty(polarAttributes)
                    for j = 1:2:size(polarAttributes,2)
                        if ~isempty(polarAttributes{i,j})
                            set(h(end),polarAttributes{i,j}, polarAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    else
        for i = 1:size(rho,1)
            h(end+1) =  polarplot(rho(i,:));
            if ~isempty(polarAttributes)
                for j = 1:2:size(polarAttributes,2)
                    if ~isempty(polarAttributes{i,j})
                        set(h(end),polarAttributes{i,j}, polarAttributes{i,j+1})
                    end
                end
            end
        end
    end
end

function h = handlePolarFittingCurves(figureStruct, or, h, mode)
global fitLineWidthDefault
if strcmpi(mode, 'normal')
    thetaFitField = 'thetaFit';
    rhoFitField = 'rhoFit';
    rhoFitLineStyle = 'polarFitLineStyle';
    rhoFitLineWidth = 'polarFitLineWidth';
    polarFitAttributes = 'polarFitAttributes';
else
    thetaFitField = 'zoomFrameThetaFit';
    rhoFitField = 'zoomFrameRhoFit';
    rhoFitLineStyle = 'zoomFitLineStyle';
    rhoFitLineWidth = 'zoomFitLineWidth';
    polarFitAttributes = 'zoomFramePolarFitAttributes';
end
if ~ isfield(figureStruct, rhoFitLineWidth)
    figureStruct.(rhoFitLineWidth) = fitLineWidthDefault;
end

if  isfield(figureStruct,thetaFitField)
    thetaFit = figureStruct.(thetaFitField);
    if or < 0
        thetaFit = thetaFit.';
    end
    
    if isfield(figureStruct, polarFitAttributes)
        fitAttributes = figureStruct.(polarFitAttributes);
        if size(fitAttributes,1)
            fitAttributes = repmat(fitAttributes,size(thetaFit,1),1);
        end
    else
        fitAttributes = {};
    end
    
    if isfield(figureStruct, rhoFitField)
        rhoFit = figureStruct.(rhoFitField);
        if or < 0
            rhoFit = rhoFit.';
        end
        if isfield(figureStruct,rhoFitLineStyle)
            lineStyle = figureStruct.(rhoFitLineStyle);
            if iscell(lineStyle)
                for i = 1:size(rhoFit,1)
                    h(end+1) = polarplot(thetaFit(i,:),rhoFit(i,:),lineStyle{i}, 'LineWidth', figureStruct.(rhoFitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(rhoFit,1)
                    h(end+1) = polarplot(thetaFit(i,:),rhoFit(i,:),lineStyle, 'LineWidth', figureStruct.(rhoFitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                h(end),fitAttributes{i,j}, fitAttributes{i,j+1}
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
            
        else
            for i = 1:size(rhoFit,1)
                h(end+1) = polarplot(thetaFit(i,:),rhoFit(i,:),'r-', 'LineWidth', figureStruct.(rhoFitLineWidth));
                if ~isempty(fitAttributes)
                    for j = 1:2:length(fitAttributes)
                        if ~isempty(fitAttributes{i,j})
                            set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    end
else
    if isfield(figureStruct, rhoFitField)
        rhoFit = figureStruct.(rhoFitField);
        if or < 0
            rhoFit = rhoFit.';
        end
        if isfield(figureStruct,rhoFitLineStyle)
            lineStyle = figureStruct.(rhoFitLineStyle);
            if iscell(lineStyle)
                for i = 1:size(rhoFit,1)
                    h(end+1) = polarplot(rhoFit(i,:),lineStyle{i}, 'LineWidth', figureStruct.(rhoFitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            else
                for i = 1:size(rhoFit,1)
                    h(end+1) = plot(rhoFit(i,:),lineStyle, 'LineWidth', figureStruct.(rhoFitLineWidth));
                    if ~isempty(fitAttributes)
                        for j = 1:2:length(fitAttributes)
                            if ~isempty(fitAttributes{i,j})
                                set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                            end
                        end
                    end
                end
            end
            
        else
            for i = 1:size(rhoFit,1)
                h(end+1) = polarplot(rhoFit(i,:),'r-', 'LineWidth', figureStruct.(rhoFitLineWidth));
                if ~isempty(fitAttributes)
                    for j = 1:2:length(fitAttributes)
                        if ~isempty(fitAttributes{i,j})
                            set(h(end),fitAttributes{i,j}, fitAttributes{i,j+1})
                        end
                    end
                end
            end
        end
    end
end

function h = handle3DPlot(figureStruct, h)

if isfield(figureStruct, 'method3D')
    method = figureStruct.method3D;
    if strcmpi(method, 'contour')
        if isfield(figureStruct, 'Z')
            Z = figureStruct.Z;
            if isfield(figureStruct, 'XY')
                XY = figureStruct.XY;
                if isfield(figureStruct, 'contourValues')
                    conVals = figureStruct.contourValues;
                    for i = 1:size(Z,3)
                        h(end+1) = contour(XY(:,:,1,i),XY(:,:,2,i),Z(:,:,i),conVals(:,i));
                    end
                else
                    for i = 1:size(Z,3)
                        h(end+1) = contour(XY(:,:,1,i),XY(:,:,2,i),Z(:,:,i));
                    end
                end
            else
                for i = 1:size(Z,3)
                    h(end+1) = contour(Z(:,:,i));
                end
            end
        end
    elseif strcmpi(method, 'surf')
        if isfield(figureStruct, 'Z')
            Z = figureStruct.Z;
            if isfield(figureStruct, 'XY')
                XY = figureStruct.XY;
                for i = 1:size(Z,3)
                    h(end+1) = surf(XY(:,:,1,i),XY(:,:,2,i),Z(:,:,i));
                end
                view(3)
            else
                for i = 1:size(Z,3)
                    h(end+1) = surf(Z(:,:,i));
                end
                view(3)
            end
        end
    end
end

function figureStruct = handleAppearance(figureStruct, mode)
if strcmp(mode,'zoom')
    ax1 = 'zoomAxes';
    xlimits = 'zoomXlimits';
    ylimits = 'zoomYlimits';
    fontS = 'zoomFontSize';
    xLabel = 'zoomXLabel';
    yLabel = 'zoomYLabel';
    zLabel = 'zoomZLabel';
    axesGrid = 'zoomAxesGrid';
    leg = 'zoomLegend';
    shouldBox = 'zoomBox';
elseif strcmp(mode,'rightaxis')
    if isempty(cell2mat(strfind(fieldnames(figureStruct), 'right')))
        return
    end
    ax1 = 'rightAxes';
    xlimits = 'rightXlimits';
    ylimits = 'rightYlimits';
    fontS = 'rightFontSize';
    xLabel = 'rightXLabel';
    yLabel = 'rightYLabel';
    zLabel = 'rightZLabel';
    axesGrid = 'rightAxesGrid';
    leg = 'rightLegend';
    shouldBox = 'rightZoomBox';
    yyaxis right
elseif strcmp(mode,'normal')
    ax1 = 'axes';
    xlimits = 'xlimits';
    ylimits = 'ylimits';
    fontS = 'fontSize';
    xLabel = 'xLabel';
    yLabel = 'yLabel';
    zLabel = 'zLabel';
    axesGrid = 'axesGrid';
    leg = 'legend';
    shouldBox = 'box';
else
    return
end
if isfield(figureStruct,fontS)
    set(gca,'FontSize',figureStruct.(fontS))
else
    if strcmp(mode, 'zoom')
        set(gca,'FontSize',10)
    else
        set(gca,'FontSize',15)
    end
end
if isfield(figureStruct,xLabel)
    xlabel(figureStruct.(xLabel),'interpreter', 'latex')
end
if isfield(figureStruct,yLabel)
    ylabel(figureStruct.(yLabel),'interpreter', 'latex')
end
if isfield(figureStruct,zLabel)
    zlabel(figureStruct.(zLabel),'interpreter', 'latex')
end
if isfield(figureStruct,axesGrid)
    if figureStruct.(axesGrid)
        grid on
    end
else
    grid off
end
if isfield(figureStruct, shouldBox)
    if figureStruct.(shouldBox)
        box on
    else
        box off
    end
else
    box on
end
if isfield(figureStruct,ax1)
    axis(figureStruct.(ax1))
end
if isfield(figureStruct, xlimits)
    xlim(figureStruct.(xlimits))
end
if isfield(figureStruct,ylimits)
    ylim(figureStruct.(ylimits))
end

if isfield(figureStruct, 'thetaLimits')
   thetalim(figureStruct.thetaLimits) 
end
if isfield(figureStruct, 'rhoLimits')
   rlim(figureStruct.rhoLimits) 
end
if isfield(figureStruct, 'thetaZeroLocation')
   ax = gca;
   ax.ThetaZeroLocation = figureStruct.thetaZeroLocation;
end

if isfield(figureStruct, leg)
    lgd = legend(figureStruct.(leg));
    if isfield(figureStruct,[leg 'Attributes'])
        att = figureStruct.([leg 'Attributes']);
        for i = 1:2:size(att, 2)
            set(lgd, att{i}, att{i+1})
        end
    end
end
if strcmpi(mode, 'rightaxis')
    yyaxis left
end


function handleSlim(figureStruct)
if isfield(figureStruct, 'slimOff') && figureStruct.slimOff
    return
end
    
if ~isfield(figureStruct, 'figuresModified')
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    if isfield(figureStruct, 'figureNumber')
        figureStruct.figuresModified = figureStruct.figureNumber;
    else
        figureStruct.figuresModified = 1;
    end
else % figuresModified exists
    if isfield(figureStruct, 'figureNumber') % der findes forskellige figurer
        if ~any(figureStruct.figuresModified == figureStruct.figureNumber)
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            figureStruct.figuresModified = [figureStruct.figuresModified figureStruct.figureNumber];
        end
    end
end

function handleZoomFrame(figureStruct, or)
if isfield(figureStruct, 'zoomFrameMethod')
    s = '';
    if strcmpi(figureStruct.zoomFrameMethod, 'Automatic')
        data = extractData(figureStruct,or);
        findBestPosition(data);
    else
        clc
        disp('You gave chosen to place the zoom frame manually')
        if exist('figurePreferences.mat','file') == 2
            while ~(strcmpi(s,'y') || strcmpi(s,'n'))
                s = input('Preference file detected - do you wish to use this position for the zoom frame? [''y''/''n'']');
            end
        end
        if strcmpi(s,'y')
            inp = load('figurePreferences.mat','preferredPosition');
            axes('Position',inp.preferredPosition)
        else
            zoomFramePositionUI();
        end
    end
    box on;
    hold on;
    [or, h] = handleErrorBars(figureStruct, 'zoom');
    h = handleFittingCurves(figureStruct, or, h, 'zoom');
    handleAppearance(figureStruct, 'zoom')
    figureStruct.zoomFramePlotHandles = h;
    figure(gcf)
end



function data = extractData(figureStruct,or)
if isfield(figureStruct, 'y')
    y = figureStruct.y;
    if or > 0
        y = y.';
    end
    data(2,:) = reshape(y,1, numel(y));
    if isfield(figureStruct, 'x')
        x = figureStruct.x;
        if or > 0
            x = x.';
        end
        data(1,:) = reshape(x,1,numel(x));
    else
        data(1,:) = 1:length(y);
    end
end
if  isfield(figureStruct,'xfit')
    xfit = figureStruct.xfit;
    if or > 0
        xfit = xfit.';
    end
    if isfield(figureStruct, 'yfit')
        or = getOrientation(figureStruct.yfit);
        yfit = figureStruct.yfit;
        if or > 0
            yfit = yfit.';
        end
        xfit = reshape(xfit, 1, numel(xfit));
        yfit = reshape(yfit, 1, numel(yfit));
        f = [xfit;yfit];
        data = [data f ];
    end
end


function moveFrame(width, height)
global pos finished
finished = 0;
pos = [0.17 0.15];
figure(gcf)
axes('Position',[pos(1) pos(2) width height])
set(gcf, 'KeyPressFcn',@handlekeyboardinput);
hold on
while ~finished
    set(gca, 'Position', [pos(1) pos(2) width height])
    pause(0.1)
end
preferredPosition = [pos(1) pos(2) width height];
tosave = '';
while ~(strcmpi(tosave,'y') || strcmpi(tosave,'n'))
    tosave = input('Do you wish to save this as your preferred position? [''y''/''n'']\n');
end
if strcmp(tosave, 'y')
    save('figurePreferences.mat', 'preferredPosition')
end

function zoomFramePositionUI()
width = -1;
height = -1;
while ~(width >=0 && width <=1)
    width = input('Please enter width of the zoom frame [0 1]\n');
end
while ~(height >=0 && height <=1)
    height =  input('Please enter height of the zoom frame [0 1]\n');
end

disp('Use "wasd" to steer the zoom frame to the desired location. Press "f" when you are finsished')
moveFrame(width, height);


function arrowsUI(figureStruct)
global finished
%disp('Click at the positions of the arrows ending with the arrow head.')
s1 = input('Do you wish to style the arrows? [''y''/''n'']');
s3 = input('Do you wish add legends the arrows? [''y''/''n'']');
s5 = input('Do you wish to add text? [''y''/''n'']');
arrows = [];
styles = {};
arrowLegends = {};
arrowHandles = [];
finished = false;
while ~finished
    %[x, y] = ginput(2);
    %dx = x(2)-x(1);
    %dy = y(2)-y(1);
    %P = [x(1) y(1) dx dy];
    %arrows = [arrows; P];
    arrowObject = arrow();
    %arrowObject = quiver(P(1),P(2),P(3),P(4));
    if strcmpi(s1,'y')
        s2success = false;
        while ~s2success
            s2 = input('Enter arrow style properties ''{''PropertyOne'',valueOne,... }'' ');
            try
                styleCell = eval(s2);
                s2success = true;
            catch
                disp('Invalid input - try again!')
            end
        end
        for j = 1:2:size(styleCell,2)
            set(arrowObject,styleCell{j}, styleCell{j+1})
            if j == 1
                indexPastEnd = 1;
            else
                indexPastEnd = 0;
            end
            styles{end+indexPastEnd,j} = styleCell{j};
            styles{end,j+1} = styleCell{j+1};
        end
    end
    if strcmpi(s3,'y')
        s4 = input('Enter legend eg. ''sin(x)'' ');
        %legendHandles(end+1) = legend(arrowObject, '-DynamicLegend', s4)
        arrowLegends{end+1} = s4;
    end
    if strcmpi(s5, 'y')
        s6 = input('Enter arrow text');
        arrowPoints = arrowObject.Vertices;
        textRelVectorX = arrowPoints(2,1) + (arrowPoints(2,1)-arrowPoints(1,1));
        textRelVectorY = arrowPoints(2,2) + (arrowPoints(2,2)-arrowPoints(1,2));
        text(textRelVectorX, textRelVectorY, s6, 'interpreter', 'latex', 'FontSize', 15)
    end
    s = input('Are you finished? [''y''/''n'']');
    if strcmpi(s,'y')
        finished = true;
    end
    arrowHandles(end+1) = arrowObject;
end
if strcmpi(s3,'y')
    disp('Legend happens!')
    legend(arrowHandles, arrowLegends{:})
end

if isfield(figureStruct, 'filename')
    savename = ['figurePreferences' figureStruct.filename '.mat'];
else
    savename = 'figurePreferences.mat';
end
save(savename, 'arrowHandles', 'styles', 'arrowLegends')


function handleArrows(figureStruct)
hold on
if isfield(figureStruct,'arrows')
    % Format is: [x1 y1 dx1 dy1; x2 y2 ...]
    arrowPoints = figureStruct.arrows;
    if isfield(figureStruct, 'arrowAttributes')
        arrowAttributes = figureStruct.arrowAttributes;
        if size(arrowAttributes,1)
            arrowAttributes = repmat(arrowAttributes,size(arrowPoints,1),1);
        end
    else
        arrowAttributes = {};
    end
    if isfield(figureStruct, 'arrowLegend')
        arrowLegend = figureStruct.arrowLegend;
    else
        arrowLegend = {};
    end
    
    for i = 1:size(arrowPoints,1)
        arrowObject = quiver(arrowPoints(i,1), arrowPoints(i,2), arrowPoints(i,3), arrowPoints(i,4));
        if ~isempty(arrowAttributes)
            for j = 1:2:size(arrowAttributes,2)
                if ~isempty(arrowAttributes{i,j})
                    set(arrowObject,arrowAttributes{i,j}, arrowAttributes{i,j+1})
                end
            end
        end
        if ~isempty(arrowLegend)
            legend(arrowObject,'-DynamicLegend', arrowLegend{i});
        end
    end
end
s = '';
if isfield(figureStruct, 'placeArrows')
    %clc
    disp('You gave chosen to place arrows manually')
    if isfield(figureStruct, 'filename')
        loadname = ['figurePreferences' figureStruct.filename '.mat'];
    else
        loadname = 'figurePreferences.mat';
    end
    if exist(loadname,'file') == 2
        
        inp = load(loadname);
        if ~isempty(fieldnames(inp))
            while ~(strcmpi(s,'y') || strcmpi(s,'n'))
                s = input('Preference file detected - do you wish to use these arrows? [''y''/''n'']');
            end
            if strcmpi(s,'y')
                arrowHandles = inp.arrowHandles;
                styles = inp.styles;
                arrowLegends = inp.arrowLegends;
            end
        end
    end
    if strcmpi(s,'y')
        hold on
        for i = 1:size(arrows,1)
            arrowObject = quiver(arrows(i,1), arrows(i,2), arrows(i,3), arrows(i,4));
            if ~isempty(styles)
                for j = 1:2:size(styles,2)
                    if ~isempty(styles{i,j})
                        set(arrowObject,styles{i,j}, styles{i,j+1})
                    end
                end
            end
            if ~isempty(arrowLegend)
                lgd = legend(arrowObject, '-DynamicLegend', arrowLegend{i});
                set(lgd,'Interpreter','latex')
            end
        end
    else
        arrowsUI(figureStruct);
    end
end

function handlekeyboardinput(~,keydata)
global pos finished
tast = keydata.Character;
if tast == 'w' || tast == 'W'
    pos = pos + [0 0.01];
elseif tast == 'a' || tast == 'A'
    pos = pos + [-0.01 0];
elseif tast == 's' || tast == 'S'
    pos = pos + [0 -0.01];
elseif tast == 'd' || tast == 'D'
    pos = pos + [0.01 0];
elseif tast == 'f' || tast == 'F'
    finished = 1;
end





