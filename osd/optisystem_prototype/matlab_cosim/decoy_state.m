%==========================================================================
% OptiSystem MATLAB Component: X-window Decoy State Generator (SNS-TF-QKD)
%--------------------------------------------------------------------------
% Purpose:
%   For every X-window, randomly selects one of three decoy states
%   (Vacuum, Decoy-1, Decoy-2) used for parameter estimation in
%   SNS-TF-QKD. Non-X-windows are passed straight through as vacuum with
%   no decoy logic performed. This block performs SELECTION and LOGGING
%   only — the actual attenuation to mu1/mu2 is applied elsewhere
%   (attenuator/MZM downstream); this component just decides which state
%   was chosen and whether the pulse is forwarded or blanked.
%
%   Decision logic per X-window:
%     r = rand;
%     r <  P0            -> VACUUM   , mu = mu0, output vacuum
%     P0 <= r < P0+P1     -> DECOY_1  , mu = mu1, output pulse unchanged
%     r >= P0+P1          -> DECOY_2  , mu = mu2, output pulse unchanged
%
% OptiSystem wiring:
%   Inputs  : Port 1 = Electrical (X-window indicator: 1 = X, 0 = not X)
%             Port 2 = Optical    (laser source pulse, pre-attenuation)
%   Outputs : Port 1 = Optical    (decoy-selected pulse, or vacuum)
%             Port 2 = Electrical (constant "1" on every window — see
%                                  below; NOT decoy logic)
%   Main tab: Sampled signal domain = Time
%   Run command: DecoyStateGenerator
%   User Parameters (all user-editable):
%     Parameter0 = P0   (vacuum probability,   default 0.2)
%     Parameter1 = P1   (decoy-1 probability,  default 0.4)
%     Parameter2 = P2   (decoy-2 probability,  default 0.4)
%     Parameter3 = mu0  (vacuum intensity,     default 0)
%     Parameter4 = mu1  (decoy-1 intensity,    default 0.05)
%     Parameter5 = mu2  (decoy-2 intensity,    default 0.20)
%
% Output Port 2 — window-count reference (not decoy-related):
%   A plain electrical "1" on every window, X or Z alike. This exists
%   solely so a downstream block with no input of its own (e.g.
%   RandomPhaseGenerator.m) can be wired to something and read off
%   numWindows = length(signal) — and, as a bonus, a correctly-timed
%   Time vector, since this port is built by copying InputPort1's own
%   Electrical structure (which already carries OptiSystem's real Time
%   base) rather than one built from scratch. Nothing about its value
%   is meaningful beyond its length; downstream blocks should not read
%   the signal content itself, only length(Signal).
%
% Modeling assumption:
%   Same one-sample-per-window convention used by the Window Selection
%   and SNS Decision blocks upstream: each element of InputPort1's
%   electrical vector is one window, and each channel of InputPort2
%   carries exactly one optical sample per window.
%
% Side effect:
%   Reads the existing alice_signal_record.csv (written by the Window
%   Selection and SNS Decision blocks) and ADDS this block's own columns
%   rather than touching anyone else's:
%     Window_Type          ("X" or "Z", stamped for every window)
%     Decoy_State           (X-windows only: VACUUM / DECOY_1 / DECOY_2)
%     Decoy_Mu               (X-windows only: mu0 / mu1 / mu2)
%     Decoy_Random_Number    (X-windows only: the rand() draw for this window)
%     P0, P1, P2              (X-windows only: logged for reproducibility)
%   The SNS Decision block's own Mu and Random_Number columns (which
%   describe its separate Z-window send decision) are never read or
%   modified by this block — the two are deliberately kept as distinct
%   columns since they represent different quantities from different
%   stages of the transmitter chain, even though both happen to be
%   "a mu" and "a random draw" in a loose sense.
%==========================================================================


outputDirectory = 'D:\drdo\osd\optisystem_prototype\matlab_cosim';
csvFileName     = fullfile(outputDirectory, 'alice_signal_record.csv');

if exist('Parameter0', 'var') && ~isempty(Parameter0), P0  = Parameter0; else, P0  = 0.2;  end
if exist('Parameter1', 'var') && ~isempty(Parameter1), P1  = Parameter1; else, P1  = 0.4;  end
if exist('Parameter2', 'var') && ~isempty(Parameter2), P2  = Parameter2; else, P2  = 0.4;  end
if exist('Parameter3', 'var') && ~isempty(Parameter3), mu0 = Parameter3; else, mu0 = 0;    end
if exist('Parameter4', 'var') && ~isempty(Parameter4), mu1 = Parameter4; else, mu1 = 0.05; end
if exist('Parameter5', 'var') && ~isempty(Parameter5), mu2 = Parameter5; else, mu2 = 0.20; end


if abs((P0 + P1 + P2) - 1) > 1e-9
    warning('DecoyStateGenerator:ProbabilitiesNotNormalized', ...
        'P0 + P1 + P2 = %.6f, expected 1.0. Check the User Parameters.', ...
        P0 + P1 + P2);
end


if ~strcmp(InputPort1.TypeSignal, 'Electrical')
    error('DecoyStateGenerator:InvalidControlSignal', ...
        'InputPort1 must be the Electrical X-window indicator signal.');
end
if ~strcmp(InputPort2.TypeSignal, 'Optical')
    error('DecoyStateGenerator:InvalidOpticalSignal', ...
        'InputPort2 must be the Optical laser-source signal.');
end

controlSignal = InputPort1.Sampled.Signal;   % 1xN: 1 = X-window, 0 = not X
numWindows    = length(controlSignal);

[ls, numChannels] = size(InputPort2.Sampled);
if ls == 0
    error('DecoyStateGenerator:EmptyOpticalSignal', ...
        'InputPort2.Sampled is empty — no optical channel to process.');
end

for c = 1:numChannels
    if size(InputPort2.Sampled(1, c).Signal, 2) ~= numWindows
        error('DecoyStateGenerator:SizeMismatch', ...
            ['Optical channel %d has %d samples but the electrical ' ...
             'control signal has %d windows. Align sample rates ' ...
             'before this block.'], c, ...
            size(InputPort2.Sampled(1, c).Signal, 2), numWindows);
    end
end


OutputPort1 = InputPort2;

OutputPort2 = InputPort1;
OutputPort2.Sampled.Signal = ones(1, numWindows);


isXWindowLog  = false(numWindows, 1);
decoyStateLog = strings(numWindows, 1);
decoyMuLog    = nan(numWindows, 1);
randomLog     = nan(numWindows, 1);


for w = 1:numWindows

    isXWindow = controlSignal(w) > 0.5;
    isXWindowLog(w) = isXWindow;

    if ~isXWindow
 
        for c = 1:numChannels
            OutputPort1.Sampled(1, c).Signal(:, w) = 0;
        end
        continue;
    end


    r = rand;
    randomLog(w) = r;

    if r < P0
        decoyStateLog(w) = "VACUUM";
        decoyMuLog(w)     = mu0;
        for c = 1:numChannels
            OutputPort1.Sampled(1, c).Signal(:, w) = 0;   % vacuum
        end

    elseif r < (P0 + P1)
        decoyStateLog(w) = "DECOY_1";
        decoyMuLog(w)     = mu1;
        for c = 1:numChannels
            OutputPort1.Sampled(1, c).Signal(:, w) = InputPort2.Sampled(1, c).Signal(:, w);
        end

    else
        decoyStateLog(w) = "DECOY_2";
        decoyMuLog(w)     = mu2;
        for c = 1:numChannels
            OutputPort1.Sampled(1, c).Signal(:, w) = InputPort2.Sampled(1, c).Signal(:, w);
        end
    end

end


if exist(csvFileName, 'file')
    logTable = readtable(csvFileName);
else
    % No prior log found from the upstream blocks — start a fresh table
    % rather than fail the simulation.
    warning('DecoyStateGenerator:MissingUpstreamLog', ...
        'alice_signal_record.csv not found at %s; creating it fresh.', ...
        csvFileName);
    logTable = table((1:numWindows)', 'VariableNames', {'Window_Index'});
end

% Ensure every column this block writes to exists, creating it if this
% is the first block in the chain to reference it.
if ~ismember('Window_Type', logTable.Properties.VariableNames)
    logTable.Window_Type = strings(height(logTable), 1);
end
if ~ismember('Decoy_State', logTable.Properties.VariableNames)
    logTable.Decoy_State = strings(height(logTable), 1);
end
if ~ismember('Decoy_Mu', logTable.Properties.VariableNames)
    logTable.Decoy_Mu = nan(height(logTable), 1);
end
if ~ismember('Decoy_Random_Number', logTable.Properties.VariableNames)
    logTable.Decoy_Random_Number = nan(height(logTable), 1);
end
if ~ismember('P0', logTable.Properties.VariableNames)
    logTable.P0 = nan(height(logTable), 1);
end
if ~ismember('P1', logTable.Properties.VariableNames)
    logTable.P1 = nan(height(logTable), 1);
end
if ~ismember('P2', logTable.Properties.VariableNames)
    logTable.P2 = nan(height(logTable), 1);
end


for w = 1:numWindows
    rowIdx = find(logTable.Window_Index == w, 1);
    if isempty(rowIdx)
        warning('DecoyStateGenerator:MissingRow', ...
            'No existing CSV row found for Window_Index %d; skipping.', w);
        continue;
    end

    if isXWindowLog(w)
        logTable.Window_Type(rowIdx)         = "X";
        logTable.Decoy_State(rowIdx)         = decoyStateLog(w);
        logTable.Decoy_Mu(rowIdx)            = decoyMuLog(w);
        logTable.Decoy_Random_Number(rowIdx) = randomLog(w);
        logTable.P0(rowIdx)                  = P0;
        logTable.P1(rowIdx)                  = P1;
        logTable.P2(rowIdx)                  = P2;
    else
        logTable.Window_Type(rowIdx) = "Z";
    end
end

writetable(logTable, csvFileName);