function PAPER_alphaGamma_OR_withDeviationHists_cleanMaps_SR450_relaxedOnly
%% ========================================================================
%  PAPER_alphaGamma_OR_withDeviationHists_cleanMaps_SR450_relaxedOnly.m
%
%  Purpose
%  -------
%  Generate a paper-style, SR450-only relaxed-cleanup output package.
%
%  This script is intended for:
%    - improved SR450 visualization
%    - SR450-only sensitivity analysis
%    - paper-style clean maps and deviation histograms
%
%  It should NOT silently replace the uniform all-sample paper pipeline.
%
%  Outputs
%  -------
%  1) phase-fraction summary CSV
%  2) austenite grain-size tables before/after deletion
%  3) austenite grain-size histograms before/after deletion
%  4) alpha/gamma OR segment table CSV
%  5) alpha/gamma OR summary CSV
%  6) normalized histograms:
%       - alpha/gamma misorientation
%       - deviation to KS
%       - deviation to NW
%       - minimum deviation
%  7) clean paper-style OR map:
%       - purple ferrite
%       - yellow austenite
%       - KS and NW only
%       - Other hidden
%       - rasterized base map to avoid white specks
%
%  Relaxed SR450 settings
%  ----------------------
%    CI               = 0.05
%    minGrainPixels   = 3
%    smoothN          = 1
%    grainTolDeg      = 5
%    OR tolerance     = 5 deg
%% ========================================================================

close all; clc;

%% ------------------------------------------------------------------------
% 0. Resolve paths and initialize MTEX
% -------------------------------------------------------------------------
thisFile = mfilename('fullpath');
if isempty(thisFile)
    error(['Save this script as ' ...
        'PAPER_alphaGamma_OR_withDeviationHists_cleanMaps_SR450_relaxedOnly.m and run again.']);
end

rootDir = fileparts(thisFile);
disp(['rootDir = ' rootDir]);

mtexDir = fullfile(rootDir,'mtex-6.1.0');
if exist(mtexDir,'dir')
    run(fullfile(mtexDir,'startup_mtex.m'));
else
    warning(['mtex-6.1.0 folder not found next to this script. ' ...
             'MTEX must already be on the MATLAB path.']);
end

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

oldDefaultFigureVisible = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

%% ------------------------------------------------------------------------
% 1. User parameters
% -------------------------------------------------------------------------
dataFile    = '03_SR450.ang';
sampleLabel = 'SR450';

% ---------- import / reference-frame ----------
pImport.edaxSetting = 'setting 2';
pImport.refConv     = 'convertEuler2SpatialReferenceFrame';
pImport.plotLikeOIM = true;

% ---------- display-only transform ----------
pDisplay.applyOIMLikeView = true;
pDisplay.flipX = false;
pDisplay.flipY = true;

% ---------- relaxed SR450 parameters ----------
pRun.ciMin          = 0.05;
pRun.grainTolDeg    = 5.0;
pRun.minGrainPixels = 3;
pRun.smoothN        = 1;
pRun.tolDeg         = 5.0;

% ---------- histogram bins ----------
pHist.misorientationEdgesDeg = 0:2:64;
pHist.deviationEdgesDeg      = 0:1:65;

% ---------- map colors ----------
pColor.ferrite    = [61 38 168] / 255;    % deep violet-purple
pColor.austenite  = [249 250 20] / 255;   % bright yellow
pColor.notIndexed = [1 1 1];              % white

fName = fullfile(rootDir,dataFile);
assert(exist(fName,'file') == 2, 'Missing file: %s', fName);

outDir = fullfile(rootDir,'PAPER_alphaGamma_OR_withDeviationHists_cleanMaps_SR450_relaxedOnly');
if exist(outDir,'dir')
    oldPng = dir(fullfile(outDir,'*.png'));
    for k = 1:numel(oldPng)
        delete(fullfile(outDir, oldPng(k).name));
    end
    oldCsv = dir(fullfile(outDir,'*.csv'));
    for k = 1:numel(oldCsv)
        delete(fullfile(outDir, oldCsv(k).name));
    end
else
    mkdir(outDir);
end

%% ------------------------------------------------------------------------
% 2. Import raw EBSD
% -------------------------------------------------------------------------
ebsd0 = EBSD.load(fName, pImport.refConv, pImport.edaxSetting);

if pImport.plotLikeOIM
    ebsd0.how2plot.east = xvector;
    ebsd0.how2plot.outOfScreen = -zvector;
end

ferriteName   = resolvePhaseName(ebsd0,'Ferrite');
austeniteName = resolvePhaseName(ebsd0,'Austenite');

if isempty(ferriteName)
    error('Ferrite phase not found in %s.', dataFile);
end
if isempty(austeniteName)
    error('Austenite phase not found in %s.', dataFile);
end

% Keep ferrite + austenite + notIndexed only
keepMask = ~ebsd0.isIndexed;
keepMask = keepMask | ismember(ebsd0.phaseId, unique(ebsd0(ferriteName).phaseId));
keepMask = keepMask | ismember(ebsd0.phaseId, unique(ebsd0(austeniteName).phaseId));
ebsdKeep = ebsd0(keepMask);

%% ------------------------------------------------------------------------
% 3. Phase-fraction summary
% -------------------------------------------------------------------------
ebsdCI = applyCIFilter(ebsdKeep, pRun.ciMin);

Tphase = table();
Tphase = [Tphase; phaseSummary(ebsd0,                 sampleLabel, "raw_import_all_points",        ferriteName, austeniteName)];
Tphase = [Tphase; phaseSummary(ebsdKeep,              sampleLabel, "keep_F_A_notIndexed",          ferriteName, austeniteName)];
Tphase = [Tphase; phaseSummary(ebsdKeep('indexed'),   sampleLabel, "indexed_only_before_CI",       ferriteName, austeniteName)];
Tphase = [Tphase; phaseSummary(ebsdCI,                sampleLabel, "after_CI_0p05",                ferriteName, austeniteName)];
Tphase = [Tphase; phaseSummary(ebsdCI('indexed'),     sampleLabel, "indexed_only_after_CI_0p05",   ferriteName, austeniteName)];

writetable(Tphase, fullfile(outDir,[sampleLabel '_phase_fraction_summary.csv']));
disp(' ');
disp('================ SR450 RELAXED PAPER PHASE SUMMARY ================');
disp(Tphase);

%% ------------------------------------------------------------------------
% 4. OR pipeline on relaxed-cleanup dataset
% -------------------------------------------------------------------------
ebsdB0 = ebsdCI('indexed');

[grainsB1, ebsdB0.grainId] = calcGrains(ebsdB0, 'angle', pRun.grainTolDeg * degree);

Ta_pre = austeniteGrainSizeTable(grainsB1, ebsdB0, austeniteName);
writetable(Ta_pre, fullfile(outDir,[sampleLabel '_austenite_grain_sizes_beforeDeletion.csv']));

TgrainPre = grainSummary(grainsB1, ebsdB0, sampleLabel, ...
    "after_first_calcGrains", ferriteName, austeniteName, pRun.minGrainPixels);
writetable(TgrainPre, fullfile(outDir,[sampleLabel '_grainSummary_beforeDeletion.csv']));
disp(' ');
disp('================ SR450 RELAXED GRAIN SUMMARY BEFORE DELETION ================');
disp(TgrainPre);

smallGrains = grainsB1(grainsB1.numPixel < pRun.minGrainPixels);
ebsdB1 = ebsdB0;
if ~isempty(smallGrains)
    ebsdB1(smallGrains) = [];
end

[grainsB2, ebsdB1.grainId] = calcGrains(ebsdB1, 'angle', pRun.grainTolDeg * degree);
grainsB2 = smooth(grainsB2, pRun.smoothN);

Ta_post = austeniteGrainSizeTable(grainsB2, ebsdB1, austeniteName);
writetable(Ta_post, fullfile(outDir,[sampleLabel '_austenite_grain_sizes_afterDeletion.csv']));

TgrainPost = grainSummary(grainsB2, ebsdB1, sampleLabel, ...
    "after_small_grain_deletion_and_recalc", ferriteName, austeniteName, pRun.minGrainPixels);
writetable(TgrainPost, fullfile(outDir,[sampleLabel '_grainSummary_afterDeletion.csv']));
disp(' ');
disp('================ SR450 RELAXED GRAIN SUMMARY AFTER DELETION ================');
disp(TgrainPost);

renderAusteniteGrainSizeHistogram(Ta_pre.NumPixel, ...
    fullfile(outDir,[sampleLabel '_austenite_grain_size_beforeDeletion']), ...
    [sampleLabel ' | austenite grain size before deletion'], ...
    pRun.minGrainPixels);

renderAusteniteGrainSizeHistogram(Ta_post.NumPixel, ...
    fullfile(outDir,[sampleLabel '_austenite_grain_size_afterDeletion']), ...
    [sampleLabel ' | austenite grain size after deletion'], ...
    pRun.minGrainPixels);

%% ------------------------------------------------------------------------
% 5. Cleaned alpha/gamma OR analysis
% -------------------------------------------------------------------------
gAF = grainsB2.boundary(austeniteName, ferriteName);

if isempty(gAF)
    warning('No ferrite-austenite boundaries found after relaxed cleanup.');
    set(0,'DefaultFigureVisible',oldDefaultFigureVisible);
    return;
end

aGrains = grainsB2(austeniteName);
fGrains = grainsB2(ferriteName);
[aSegID, fSegID] = resolveBoundaryPhaseIDs(gAF, aGrains, fGrains);

csA = ebsdB1(austeniteName).CS;
csF = ebsdB1(ferriteName).CS;

KS = orientation.KurdjumovSachs(csA, csF);
NW = orientation.NishiyamaWassermann(csA, csF);

mori   = gAF.misorientation;
segLen = segLength(gAF);

thetaAG = angle(mori) ./ degree;
devKS   = angle(mori, KS) ./ degree;
devNW   = angle(mori, NW) ./ degree;
devMin  = min(devKS, devNW);

[orClass, isKS, isNW, isOther] = classifyORNearest(devKS, devNW, pRun.tolDeg);

[x1, y1, x2, y2] = extractBoundarySegmentsTrue(gAF);

Tseg = table( ...
    repmat(string(sampleLabel), length(gAF), 1), ...
    (1:length(gAF))', ...
    aSegID(:), fSegID(:), ...
    thetaAG(:), segLen(:), ...
    devKS(:), devNW(:), devMin(:), orClass(:), ...
    x1(:), y1(:), x2(:), y2(:), ...
    'VariableNames', {'Sample','SegmentID','AGrainID','FGrainID', ...
    'AlphaGammaMisorientation_deg','Length_um', ...
    'DevKS_deg','DevNW_deg','DevMin_deg','ORClass', ...
    'x1_um','y1_um','x2_um','y2_um'});

writetable(Tseg, fullfile(outDir,[sampleLabel '_alphaGamma_segment_table.csv']));

Tsum = summarizeAlphaGammaOR(sampleLabel, thetaAG, segLen, isKS, isNW, isOther, devKS, devNW, devMin);
writetable(Tsum, fullfile(outDir,[sampleLabel '_alphaGamma_OR_summary.csv']));

disp(' ');
disp('================ SR450 RELAXED OR SUMMARY ================');
disp(Tsum);

%% ------------------------------------------------------------------------
% 6. Histogram CSVs and plots
% -------------------------------------------------------------------------
writeHistogramCSV(thetaAG, pHist.misorientationEdgesDeg, ...
    fullfile(outDir,[sampleLabel '_alphaGamma_misorientation_histogram.csv']));
writeHistogramCSV(devKS, pHist.deviationEdgesDeg, ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devKS_histogram.csv']));
writeHistogramCSV(devNW, pHist.deviationEdgesDeg, ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devNW_histogram.csv']));
writeHistogramCSV(devMin, pHist.deviationEdgesDeg, ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devMin_histogram.csv']));

renderHistogramFromCSV( ...
    fullfile(outDir,[sampleLabel '_alphaGamma_misorientation_histogram.csv']), ...
    fullfile(outDir,[sampleLabel '_alphaGamma_misorientation_histogram']), ...
    [sampleLabel ' | \alpha/\gamma misorientation histogram'], ...
    '\alpha/\gamma misorientation angle (deg)', ...
    0.45, ...
    [0 0.1 0.2 0.3 0.4], ...
    {'0','0.1','0.2','0.3','0.4'});

renderHistogramFromCSV( ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devKS_histogram.csv']), ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devKS_histogram']), ...
    [sampleLabel ' | deviation from K-S'], ...
    'Deviation from K-S (deg)', ...
    0.25, ...
    [0 0.05 0.10 0.15 0.20 0.25], ...
    {'0','0.05','0.10','0.15','0.20','0.25'});

renderHistogramFromCSV( ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devNW_histogram.csv']), ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devNW_histogram']), ...
    [sampleLabel ' | deviation from N-W'], ...
    'Deviation from N-W (deg)', ...
    0.25, ...
    [0 0.05 0.10 0.15 0.20 0.25], ...
    {'0','0.05','0.10','0.15','0.20','0.25'});

renderHistogramFromCSV( ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devMin_histogram.csv']), ...
    fullfile(outDir,[sampleLabel '_alphaGamma_devMin_histogram']), ...
    [sampleLabel ' | minimum deviation to nearest rational OR'], ...
    'Minimum deviation, min(\Delta\theta_{KS}, \Delta\theta_{NW}) (deg)', ...
    0.3, ...
    [0 0.1 0.2 0.3], ...
    {'0','0.1','0.2','0.3'});

%% ------------------------------------------------------------------------
% 7. Clean paper-style map: KS/NW only, Other hidden
% -------------------------------------------------------------------------
renderCleanORMapTrue(ebsdB1, x1, y1, x2, y2, isKS, isNW, isOther, ...
    ferriteName, austeniteName, outDir, ...
    [sampleLabel '_alphaGamma_OR_map_trueBoundary_clean_relaxedOnly'], ...
    [sampleLabel ' | relaxed \alpha/\gamma OR map (KS/NW only)'], ...
    pDisplay, pColor);

%% ------------------------------------------------------------------------
% 8. Key summary table
% -------------------------------------------------------------------------
Tkey = table( ...
    string(sampleLabel), ...
    Tphase.AusteniteFrac_ofIndexed(Tphase.Stage=="indexed_only_before_CI"), ...
    Tphase.AusteniteFrac_ofIndexed(Tphase.Stage=="indexed_only_after_CI_0p05"), ...
    TgrainPre.ThresholdPx, ...
    TgrainPre.nAusteniteGrains_ltThresholdPx, ...
    TgrainPost.AusteniteFrac_ofIndexed, ...
    Tsum.nAlphaGammaSegments, ...
    'VariableNames', {'Sample', ...
    'GammaFrac_beforeCI', ...
    'GammaFrac_afterCI_0p05', ...
    'ThresholdPx', ...
    'nAusteniteGrains_ltThresholdPx_beforeDeletion', ...
    'GammaFrac_afterDeletion', ...
    'nAlphaGammaSegments'});

writetable(Tkey, fullfile(outDir,[sampleLabel '_key_summary_table.csv']));
disp(' ');
disp('================ SR450 RELAXED KEY SUMMARY ================');
disp(Tkey);

disp(' ');
disp(['All outputs are in: ' outDir]);

set(0,'DefaultFigureVisible',oldDefaultFigureVisible);

end

%% ========================================================================
% Helper functions
%% ========================================================================

function ebsdOut = applyCIFilter(ebsdIn, ciMin)

ebsdOut = ebsdIn;
if isfield(ebsdOut.prop,'ci')
    lowCI = ebsdOut.isIndexed & (ebsdOut.prop.ci < ciMin);
    ebsdOut(lowCI) = [];
else
    warning('No CI field found. CI filtering not applied.');
end

end

function phaseName = resolvePhaseName(ebsd, phaseQuery)

names = cellstr(ebsd.mineralList(:));
names = names(~cellfun(@isempty,names));
phaseName = '';

if isempty(names)
    return;
end

idx = find(strcmpi(names, phaseQuery), 1, 'first');
if isempty(idx)
    idx = find(contains(lower(string(names)), lower(string(phaseQuery))), 1, 'first');
end

if ~isempty(idx)
    phaseName = names{idx};
end

end

function T = phaseSummary(ebsd, sampleLabel, stageName, ferriteName, austeniteName)

nTotal     = double(length(ebsd));
nIndexed   = double(length(ebsd('indexed')));
nFerrite   = double(length(ebsd(ferriteName)));
nAustenite = double(length(ebsd(austeniteName)));
nNotIdx    = nTotal - nIndexed;

if nIndexed > 0
    ferriteFracIndexed   = nFerrite   / nIndexed;
    austeniteFracIndexed = nAustenite / nIndexed;
else
    ferriteFracIndexed   = NaN;
    austeniteFracIndexed = NaN;
end

if nTotal > 0
    ferriteFracTotal   = nFerrite   / nTotal;
    austeniteFracTotal = nAustenite / nTotal;
    notIdxFracTotal    = nNotIdx    / nTotal;
else
    ferriteFracTotal   = NaN;
    austeniteFracTotal = NaN;
    notIdxFracTotal    = NaN;
end

T = table( ...
    string(sampleLabel), string(stageName), ...
    nTotal, nIndexed, nNotIdx, nFerrite, nAustenite, ...
    ferriteFracIndexed, austeniteFracIndexed, ...
    ferriteFracTotal, austeniteFracTotal, notIdxFracTotal, ...
    'VariableNames', {'Sample','Stage', ...
    'nTotalPoints','nIndexedPoints','nNotIndexedPoints','nFerritePoints','nAustenitePoints', ...
    'FerriteFrac_ofIndexed','AusteniteFrac_ofIndexed', ...
    'FerriteFrac_ofTotal','AusteniteFrac_ofTotal','NotIndexedFrac_ofTotal'});

end

function T = grainSummary(grains, ebsd, sampleLabel, stageName, ferriteName, austeniteName, thresholdPx)

gFerrite   = grains(ferriteName);
gAustenite = grains(austeniteName);

nGrainsTotal = double(length(grains));
nFerriteGr   = double(length(gFerrite));
nAusteniteGr = double(length(gAustenite));

nTotal     = double(length(ebsd));
nIndexed   = double(length(ebsd('indexed')));
nFerrite   = double(length(ebsd(ferriteName)));
nAustenite = double(length(ebsd(austeniteName)));
nNotIdx    = nTotal - nIndexed;

if nIndexed > 0
    ferriteFracIndexed   = nFerrite   / nIndexed;
    austeniteFracIndexed = nAustenite / nIndexed;
else
    ferriteFracIndexed   = NaN;
    austeniteFracIndexed = NaN;
end

if nTotal > 0
    ferriteFracTotal   = nFerrite   / nTotal;
    austeniteFracTotal = nAustenite / nTotal;
    notIdxFracTotal    = nNotIdx    / nTotal;
else
    ferriteFracTotal   = NaN;
    austeniteFracTotal = NaN;
    notIdxFracTotal    = NaN;
end

if ~isempty(gAustenite)
    medianAGrainPx = double(median(gAustenite.numPixel));
    meanAGrainPx   = double(mean(gAustenite.numPixel));
    nA_small       = double(sum(gAustenite.numPixel < thresholdPx));
else
    medianAGrainPx = NaN;
    meanAGrainPx   = NaN;
    nA_small       = 0;
end

T = table( ...
    string(sampleLabel), string(stageName), ...
    nTotal, nIndexed, nNotIdx, nFerrite, nAustenite, ...
    ferriteFracIndexed, austeniteFracIndexed, ...
    ferriteFracTotal, austeniteFracTotal, notIdxFracTotal, ...
    nGrainsTotal, nFerriteGr, nAusteniteGr, ...
    medianAGrainPx, meanAGrainPx, thresholdPx, nA_small, ...
    'VariableNames', {'Sample','Stage', ...
    'nTotalPoints','nIndexedPoints','nNotIndexedPoints','nFerritePoints','nAustenitePoints', ...
    'FerriteFrac_ofIndexed','AusteniteFrac_ofIndexed', ...
    'FerriteFrac_ofTotal','AusteniteFrac_ofTotal','NotIndexedFrac_ofTotal', ...
    'nGrainsTotal','nFerriteGrains','nAusteniteGrains', ...
    'MedianAusteniteGrainPx','MeanAusteniteGrainPx', ...
    'ThresholdPx','nAusteniteGrains_ltThresholdPx'});

end

function T = austeniteGrainSizeTable(grains, ebsdIndexed, austeniteName)

gA = grains(austeniteName);

if isempty(gA)
    T = table([], [], [], 'VariableNames', {'GrainID','NumPixel','Area_um2'});
    return;
end

pixelArea = estimatePixelArea(ebsdIndexed);
area_um2 = double(gA.numPixel(:)) * pixelArea;

T = table( ...
    double(gA.id(:)), ...
    double(gA.numPixel(:)), ...
    area_um2(:), ...
    'VariableNames', {'GrainID','NumPixel','Area_um2'});

end

function pixelArea = estimatePixelArea(ebsdIndexed)

x = double(ebsdIndexed.x(:));
y = double(ebsdIndexed.y(:));

ux = unique(sort(x));
uy = unique(sort(y));

if numel(ux) > 1
    dx = median(diff(ux));
else
    dx = 1.0;
end

if numel(uy) > 1
    dy = median(diff(uy));
else
    dy = 1.0;
end

if ~isfinite(dx) || dx <= 0
    dx = 1.0;
end
if ~isfinite(dy) || dy <= 0
    dy = 1.0;
end

pixelArea = dx * dy;

end

function [aSegID, fSegID] = resolveBoundaryPhaseIDs(gAF, aGrains, fGrains)

g1 = gAF.grainId(:,1);
g2 = gAF.grainId(:,2);

aIDs = aGrains.id(:);
fIDs = fGrains.id(:);

isCol1A = ismember(g1, aIDs) & ismember(g2, fIDs);
isCol2A = ismember(g2, aIDs) & ismember(g1, fIDs);

if ~all(isCol1A | isCol2A)
    error('Could not resolve austenite/ferrite grain IDs on alpha/gamma boundaries.');
end

aSegID = zeros(length(gAF),1);
fSegID = zeros(length(gAF),1);

aSegID(isCol1A) = g1(isCol1A);
fSegID(isCol1A) = g2(isCol1A);

aSegID(isCol2A) = g2(isCol2A);
fSegID(isCol2A) = g1(isCol2A);

end

function [orClass, isKS, isNW, isOther] = classifyORNearest(devKS, devNW, tolDeg)

isKS = (devKS <= tolDeg) & (devKS < devNW);
isNW = (devNW <= tolDeg) & (devNW < devKS);

tieMask = (devKS <= tolDeg) & (devNW <= tolDeg) & (abs(devKS - devNW) < 1e-12);
isKS(tieMask) = true;

isOther = ~(isKS | isNW);

orClass = strings(size(devKS));
orClass(isKS)    = "KS";
orClass(isNW)    = "NW";
orClass(isOther) = "Other";

end

function [x1, y1, x2, y2] = extractBoundarySegmentsTrue(gB)

nSeg = length(gB);

if nSeg == 0
    x1 = zeros(0,1); y1 = zeros(0,1);
    x2 = zeros(0,1); y2 = zeros(0,1);
    return;
end

mp = gB.midPoint;
mx = double(mp.x(:));
my = double(mp.y(:));

dir = gB.direction;
dx = double(dir.x(:));
dy = double(dir.y(:));

L = double(segLength(gB));
L = L(:);

dnorm = sqrt(dx.^2 + dy.^2);
dx = dx ./ dnorm;
dy = dy ./ dnorm;

halfL = 0.5 .* L;

x1 = mx - halfL .* dx;
y1 = my - halfL .* dy;
x2 = mx + halfL .* dx;
y2 = my + halfL .* dy;

end

function T = summarizeAlphaGammaOR(sampleLabel, thetaAG, segLen, isKS, isNW, isOther, devKS, devNW, devMin)

nSeg = numel(thetaAG);
totLen = sum(segLen);

nKS = sum(isKS);
nNW = sum(isNW);
nOther = sum(isOther);

lenKS = sum(segLen(isKS));
lenNW = sum(segLen(isNW));
lenOther = sum(segLen(isOther));

T = table( ...
    string(sampleLabel), ...
    double(nSeg), double(totLen), ...
    double(median(thetaAG)), ...
    double(nKS), double(nNW), double(nOther), ...
    double(nKS/nSeg), double(nNW/nSeg), double(nOther/nSeg), ...
    double(lenKS/totLen), double(lenNW/totLen), double(lenOther/totLen), ...
    double(median(devKS)), double(median(devNW)), double(median(devMin)), ...
    'VariableNames', {'Sample','nAlphaGammaSegments','TotalAlphaGammaLen_um', ...
    'MedianAlphaGammaAngle_deg', ...
    'nKS','nNW','nOther', ...
    'CountFracKS','CountFracNW','CountFracOther', ...
    'LenFracKS','LenFracNW','LenFracOther', ...
    'MedianDevKS_deg','MedianDevNW_deg','MedianDevMin_deg'});

end

function writeHistogramCSV(valuesDeg, edgesDeg, outCsv)

counts = histcounts(valuesDeg, edgesDeg, 'Normalization', 'probability');

binLeft = edgesDeg(1:end-1).';
binRight = edgesDeg(2:end).';
binCenter = (binLeft + binRight) / 2;

Th = table(binLeft, binRight, binCenter, counts(:), ...
    'VariableNames', {'BinLeft_deg','BinRight_deg','BinCenter_deg','NumberFraction'});

writetable(Th, outCsv);

end
%=================================================================================================================
function renderHistogramFromCSV(inCsv, outStem, ttl, xlab, yMaxForced, yTickValsForced, yTickLabelsForced)

if nargin < 5 || isempty(yMaxForced)
    yMaxForced = [];
end
if nargin < 6
    yTickValsForced = [];
end
if nargin < 7
    yTickLabelsForced = {};
end

if ~exist(inCsv,'file')
    return;
end

T = readtable(inCsv);
if isempty(T)
    return;
end

if isempty(yMaxForced)
    globalYMax = max(T.NumberFraction);
    globalYMax = ceil((1.05 * globalYMax) / 0.05) * 0.05;
    if globalYMax <= 0
        globalYMax = 0.10;
    end
else
    globalYMax = yMaxForced;
end

if isempty(yTickValsForced)
    if globalYMax <= 0.30
        yStep = 0.05;
    else
        yStep = 0.10;
    end
    yTickVals = 0:yStep:globalYMax;
else
    yTickVals = yTickValsForced;
end

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 900 700]);
ax = axes('Parent',fig, 'Position',[0.13 0.14 0.82 0.78]);
hold(ax,'on');

bar(ax, T.BinCenter_deg, T.NumberFraction, 1.0, ...
    'FaceColor',[0.20 0.35 0.80], ...
    'EdgeColor','k', ...
    'LineWidth',1.8);

xlabel(ax, xlab, 'FontWeight','bold', 'FontSize',22, 'Interpreter','tex');
ylabel(ax, 'Number fraction', 'FontWeight','bold', 'FontSize',22);
title(ax, ttl, 'FontWeight','bold', 'FontSize',22, 'Interpreter','tex');

set(ax, ...
    'FontWeight','bold', ...
    'FontSize',18, ...
    'LineWidth',2.5, ...
    'Box','on', ...
    'Layer','top', ...
    'TickDir','in');

xlim(ax, [min(T.BinLeft_deg) max(T.BinRight_deg)]);
ylim(ax, [0 globalYMax]);
yticks(ax, yTickVals);

if ~isempty(yTickLabelsForced)
    yticklabels(ax, yTickLabelsForced);
else
    yticklabels(ax, makePrettyTickLabels(yTickVals));
end

grid(ax, 'on');
ax.GridAlpha = 0.25;

disp(['Saved histogram: ' outStem]);
disp('YTick values used:');
disp(get(ax,'YTick'));
disp('YTick labels used:');
disp(get(ax,'YTickLabel'));

saveFigureTrue(fig, outStem);
close(fig);

end
%==========================================================================================================
function labels = makePrettyTickLabels(vals)

labels = cell(size(vals));

for i = 1:numel(vals)
    v = vals(i);

    if abs(v) < 1e-12
        labels{i} = '0';
    elseif abs(v*10 - round(v*10)) < 1e-12
        labels{i} = sprintf('%.1f', v);
    else
        labels{i} = sprintf('%.2f', v);
    end
end

end

function renderAusteniteGrainSizeHistogram(numPix, outStem, ttl, thresholdPx)

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 760 560]);
ax = axes('Parent',fig, 'Position',[0.12 0.14 0.82 0.76]);
hold(ax,'on');

if isempty(numPix)
    text(ax, 0.5, 0.5, 'No austenite grains detected', ...
        'Units','normalized', 'HorizontalAlignment','center', ...
        'FontWeight','bold');
    axis(ax,'off');
else
    edges = 0.5:1:max([20; numPix(:)+1]);
    histogram(ax, numPix, edges, ...
        'Normalization','probability', ...
        'FaceColor',[0.20 0.35 0.80], ...
        'EdgeColor','k', ...
        'LineWidth',1.0);

    xline(ax, thresholdPx, '--r', sprintf('%d px threshold', thresholdPx), ...
        'LineWidth',1.3, 'FontWeight','bold', 'LabelVerticalAlignment','middle');

    xlabel(ax,'Austenite grain size (pixels)', 'FontWeight','bold');
    ylabel(ax,'Number fraction', 'FontWeight','bold');
    set(ax,'FontWeight','bold', 'FontSize',11, 'LineWidth',1.2, 'Box','on');
    grid(ax,'on');
end

title(ax, ttl, 'FontWeight','bold', 'FontSize',13, 'Interpreter','none');
saveFigureTrue(fig, outStem);
close(fig);

end

function renderPhaseMap(ebsd, ferriteName, austeniteName, outStem, ttl, pDisplay, pColor)

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 900 760]);
ax = axes('Parent',fig, 'Position',[0.08 0.08 0.84 0.84]);
hold(ax,'on');

if any(~ebsd.isIndexed)
    en = ebsd(~ebsd.isIndexed);
    [xn, yn] = transformDisplayCoords(en.x, en.y, ebsd, pDisplay);
    scatter(ax, xn, yn, 14, repmat(pColor.notIndexed, length(en), 1), ...
        's', 'filled', 'MarkerEdgeColor','none');
end

ef = ebsd(ferriteName);
if ~isempty(ef)
    [xf, yf] = transformDisplayCoords(ef.x, ef.y, ebsd, pDisplay);
    scatter(ax, xf, yf, 14, repmat(pColor.ferrite, length(ef), 1), ...
        's', 'filled', 'MarkerEdgeColor','none');
end

ea = ebsd(austeniteName);
if ~isempty(ea)
    [xa, ya] = transformDisplayCoords(ea.x, ea.y, ebsd, pDisplay);
    scatter(ax, xa, ya, 14, repmat(pColor.austenite, length(ea), 1), ...
        's', 'filled', 'MarkerEdgeColor','none');
end

[xAll, yAll] = transformDisplayCoords(ebsd.x, ebsd.y, ebsd, pDisplay);
xlim(ax,[min(xAll) max(xAll)]);
ylim(ax,[min(yAll) max(yAll)]);
axis(ax,'equal');
set(ax,'YDir','normal');

ax.Visible = 'off';
title(ax, ttl, 'FontWeight','bold', 'FontSize',14, 'Interpreter','none');

saveFigureTrue(fig, outStem);
close(fig);

end

function renderCleanORMapTrue(ebsd, x1, y1, x2, y2, isKS, isNW, isOther, ...
    ferriteName, austeniteName, outDir, fileStem, ttl, pDisplay, pColor) %#ok<INUSD>

[fig, ax] = createSidebarFigure();

[imgRGB, xlimData, ylimData] = buildRasterPhaseMapRGB( ...
    ebsd, ferriteName, austeniteName, pDisplay, pColor);

image(ax, 'XData', xlimData, 'YData', ylimData, 'CData', imgRGB);
set(ax, 'YDir', 'normal');
hold(ax, 'on');

colKS = [0.85 0.33 0.10];
colNW = [0.00 0.45 0.74];

[x1p, y1p] = transformDisplayCoords(x1, y1, ebsd, pDisplay);
[x2p, y2p] = transformDisplayCoords(x2, y2, ebsd, pDisplay);

drawSegments(ax, x1p, y1p, x2p, y2p, isNW, colNW, 1.20);
drawSegments(ax, x1p, y1p, x2p, y2p, isKS, colKS, 1.20);

xlim(ax, xlimData);
ylim(ax, ylimData);
axis(ax, 'equal');
set(ax, 'YDir', 'normal');

ax.Visible = 'off';
ax.LineWidth = 1.2;
set(ax, 'Color', pColor.ferrite);

title(ax, ttl, 'FontWeight','bold', 'FontSize',14, 'Interpreter','tex');

addORSidebarClean(fig, colKS, colNW);
addScaleBarSidebar(fig, 40);

saveFigureTrue(fig, fullfile(outDir,fileStem));
close(fig);

end

function [imgRGB, xlimData, ylimData] = buildRasterPhaseMapRGB( ...
    ebsd, ferriteName, austeniteName, pDisplay, pColor)

eIdx = ebsd('indexed');
if isempty(eIdx)
    error('No indexed EBSD points available for raster phase map rendering.');
end

[xp, yp] = transformDisplayCoords(eIdx.x, eIdx.y, ebsd, pDisplay);

xmin = min(xp);
xmax = max(xp);
ymin = min(yp);
ymax = max(yp);

ux = unique(sort(xp));
uy = unique(sort(yp));

if numel(ux) > 1
    dx = median(diff(ux));
else
    dx = 1;
end

if numel(uy) > 1
    dy = median(diff(uy));
else
    dy = 1;
end

if ~isfinite(dx) || dx <= 0
    dx = 1;
end
if ~isfinite(dy) || dy <= 0
    dy = 1;
end

gx = xmin:dx:xmax;
gy = ymin:dy:ymax;

Nx = numel(gx);
Ny = numel(gy);

aPhaseIds = unique(ebsd(austeniteName).phaseId);
isA = ismember(eIdx.phaseId, aPhaseIds);
phaseVal = double(isA(:));

XY = [xp(:), yp(:)];
[XYu, ~, ic] = unique(XY, 'rows', 'stable');
phaseValU = accumarray(ic, phaseVal, [], @mean);
phaseValU = double(phaseValU > 0.5);

F = scatteredInterpolant(XYu(:,1), XYu(:,2), phaseValU, 'nearest', 'nearest');
[Xg, Yg] = meshgrid(gx, gy);
phaseGrid = F(Xg, Yg);

isAgrid = phaseGrid > 0.5;

imgRGB = zeros(Ny, Nx, 3);
for c = 1:3
    chan = pColor.ferrite(c) * ones(Ny, Nx);
    chan(isAgrid) = pColor.austenite(c);
    imgRGB(:,:,c) = chan;
end

xlimData = [gx(1) - dx/2, gx(end) + dx/2];
ylimData = [gy(1) - dy/2, gy(end) + dy/2];

end

function [fig, ax] = createSidebarFigure()

fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position',[100 100 1200 760]);
ax = axes('Parent',fig, 'Position',[0.06 0.08 0.66 0.84]);
hold(ax,'on');

end

function addORSidebarClean(fig, colKS, colNW)

x0 = 0.78;
y0 = 0.82;
dy = 0.065;
linew = 0.030;

annotation(fig,'textbox',[0.76 0.88 0.18 0.04], ...
    'String','\alpha/\gamma OR class', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'Interpreter','tex');

annotation(fig,'line',[x0 x0+linew],[y0 y0], ...
    'Color',colKS, 'LineWidth',2.8);
annotation(fig,'textbox',[x0+0.04 y0-0.015 0.12 0.03], ...
    'String','KS', ...
    'LineStyle','none', 'FontSize',10, 'FontWeight','bold');

annotation(fig,'line',[x0 x0+linew],[y0-dy y0-dy], ...
    'Color',colNW, 'LineWidth',2.8);
annotation(fig,'textbox',[x0+0.04 y0-dy-0.015 0.12 0.03], ...
    'String','NW', ...
    'LineStyle','none', 'FontSize',10, 'FontWeight','bold');

annotation(fig,'textbox',[0.76 0.56 0.18 0.04], ...
    'String','Base map = phase colors', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',10);

end

function addScaleBarSidebar(fig, L)

annotation(fig,'textbox',[0.76 0.28 0.16 0.04], ...
    'String','Scale bar', ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',12);

annotation(fig,'line',[0.78 0.88],[0.24 0.24], ...
    'Color','k', 'LineWidth',3.0);
annotation(fig,'textbox',[0.79 0.19 0.12 0.03], ...
    'String',sprintf('%g \\mum', L), ...
    'LineStyle','none', ...
    'FontWeight','bold', ...
    'FontSize',11);

end

function [xp, yp] = transformDisplayCoords(x, y, ebsd, pDisplay)

xp = double(x(:));
yp = double(y(:));

xmin = min(ebsd.x);
xmax = max(ebsd.x);
ymin = min(ebsd.y);
ymax = max(ebsd.y);

if pDisplay.applyOIMLikeView
    if pDisplay.flipX
        xp = xmax - (xp - xmin);
    end
    if pDisplay.flipY
        yp = ymax - (yp - ymin);
    end
end

end

function drawSegments(ax, x1, y1, x2, y2, mask, colorRGB, lineWidth)

idx = find(mask);
if isempty(idx)
    return;
end

X = nan(3*numel(idx),1);
Y = nan(3*numel(idx),1);

X(1:3:end) = x1(idx);
X(2:3:end) = x2(idx);

Y(1:3:end) = y1(idx);
Y(2:3:end) = y2(idx);

line(ax, X, Y, 'Color', colorRGB, 'LineWidth', lineWidth);

end

function saveFigureTrue(fig, fileStem)

drawnow;
pause(0.2);

try
    exportgraphics(fig, [fileStem '.png'], 'Resolution', 500);
    return;
catch
end

try
    set(fig,'PaperPositionMode','auto');
    print(fig, [fileStem '.png'], '-dpng', '-r500');
    return;
catch
end

warning('Figure export failed for: %s', fileStem);

end

