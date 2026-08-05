%==========================================================================
% OptiSystem MATLAB Component: Random Phase Generator for SNS-TF-QKD
%--------------------------------------------------------------------------
% Purpose:
%   Generates one uniformly distributed random phase phi in [0, 2*pi) per
%   simulation window and outputs it as an electrical control signal
%   intended to drive a downstream OptiSystem Phase Modulator. Also
%   assigns each phase to one of M discrete phase slices and logs the
%   result.
%
%   SCOPE NOTE: the request that named this component called it an
%   "X-window" phase generator, but the behavior actually specified (and
%   implemented here) generates a phase for EVERY window, with no
%   distinction between X and Z windows — and this component has no
%   input port carrying window-type info, so it has no way to know a
%   window's type even if it needed to. If phase randomization should
%   really be restricted to X-windows only (e.g. holding phi = 0 on
%   Z-windows), this needs an Electrical input port carrying the
%   window-type indicator, the same way AliceWindowRouter.m and
%   SNSDecisionBlock.m consume it.
%
%   CRITICAL: exactly one rand() draw happens per WINDOW, not per sample.
%   The output vector has one element per window (the same
%   one-sample-per-window convention used by every other block in this
%   pipeline), so phase stays constant across an entire pulse and will
%   not fluctuate mid-pulse and destroy TF-QKD interference.
%
% Inputs:
%   Input Port 1: Electrical, window-count/Time-base reference only —
%                 wire this to DecoyStateGenerator.m's Output Port 2
%                 (a constant "1" on every window). Only length(Signal)
%                 and the Time vector are used; the signal's actual
%                 values are ignored. The window spacing used to build
%                 this component's own output Time vector is derived
%                 directly from InputPort1.Sampled.Time (the real
%                 simulation clock) — no manually-set WindowDuration
%                 parameter is needed or used.
%
% Outputs:
%   Output Port 1: Electrical signal, phi in radians, one sample per
%                  window, wired directly into the Phase Modulator.
%
% OptiSystem wiring:
%   Inputs  : Port 1 = Electrical (window-count reference; see above)
%   Outputs : Port 1 = Electrical
%   Main tab: Sampled signal domain = Time
%   Run command: RandomPhaseGenerator
%   User Parameters (all user-editable):
%     Parameter0 = M   (number of phase slices, default 16)
%
% Side effect:
%   Reads (or creates) the shared alice_signal_record.csv and adds this
%   block's own columns, matched by Window_Index like every other block
%   in the pipeline: Phase_Rad, Phase_Deg, Slice_Number, M.
%==========================================================================

outputDirectory = 'D:\drdo\osd\optisystem_prototype\outputs';
csvFileName     = fullfile(outputDirectory, 'alice_signal_record.csv');

if exist('Parameter0', 'var') && ~isempty(Parameter0), M = Parameter0; else, M = 16; end

M = round(M);
if M < 1
    error('RandomPhaseGenerator:InvalidM', 'M must be a positive integer, got %g.', M);
end


if ~strcmp(InputPort1.TypeSignal, 'Electrical')
    error('RandomPhaseGenerator:InvalidReferenceSignal', ...
        'InputPort1 must be the Electrical window-count reference signal.');
end

refTime    = InputPort1.Sampled.Time;
numWindows = length(InputPort1.Sampled.Signal);
if numWindows < 1
    error('RandomPhaseGenerator:EmptyReferenceSignal', ...
        'InputPort1 carries zero samples — nothing to generate phases for.');
end

DeltaPhi = 2 * pi / M;   


phiLog   = zeros(1, numWindows);
sliceLog = zeros(1, numWindows);

for w = 1:numWindows

    phi = 2 * pi * rand;  

    sliceNumber = floor(phi / DeltaPhi) + 1;

    sliceNumber = min(sliceNumber, M);

    phiLog(w)   = phi;
    sliceLog(w) = sliceNumber;

end


timeVector = refTime;

OutputPort1.TypeSignal      = 'Electrical';
OutputPort1.Sampled.Signal  = phiLog;      % 1 x numWindows, radians
OutputPort1.Sampled.Time    = timeVector;


OutputPort1.Noise.Signal    = zeros(1, numWindows);
OutputPort1.Noise.Time      = timeVector;

OutputPort1.IndividualSample = [];

if exist(csvFileName, 'file')
    logTable = readtable(csvFileName);
else
    warning('RandomPhaseGenerator:MissingUpstreamLog', ...
        'alice_signal_record.csv not found at %s; creating it fresh.', ...
        csvFileName);
    logTable = table((1:numWindows)', 'VariableNames', {'Window_Index'});
end

if ~ismember('Phase_Rad', logTable.Properties.VariableNames)
    logTable.Phase_Rad = nan(height(logTable), 1);
end
if ~ismember('Phase_Deg', logTable.Properties.VariableNames)
    logTable.Phase_Deg = nan(height(logTable), 1);
end
if ~ismember('Slice_Number', logTable.Properties.VariableNames)
    logTable.Slice_Number = nan(height(logTable), 1);
end
if ~ismember('M', logTable.Properties.VariableNames)
    logTable.M = nan(height(logTable), 1);
end

for w = 1:numWindows
    rowIdx = find(logTable.Window_Index == w, 1);
    if isempty(rowIdx)
        warning('RandomPhaseGenerator:MissingRow', ...
            'No existing CSV row found for Window_Index %d; skipping.', w);
        continue;
    end

    logTable.Phase_Rad(rowIdx)    = phiLog(w);
    logTable.Phase_Deg(rowIdx)    = phiLog(w) * (180 / pi);
    logTable.Slice_Number(rowIdx) = sliceLog(w);
    logTable.M(rowIdx)            = M;
end

writetable(logTable, csvFileName);