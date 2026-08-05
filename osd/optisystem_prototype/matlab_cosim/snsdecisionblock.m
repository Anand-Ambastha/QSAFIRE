%==========================================================================
% OptiSystem MATLAB Component: SNS (Sending-or-Not-Sending) Decision Block
%--------------------------------------------------------------------------
% Purpose:
%   Implements the SNS decision process for SNS-TF-QKD. Basis (Z/X window)
%   selection has already been made upstream by the Window Selection
%   block; this component only decides, per Z-window, whether the
%   prepared coherent pulse is physically sent or withheld (vacuum sent
%   instead). X-windows are always blocked without randomization.
%
%   Truth table:
%     Electrical Input | Condition  | Optical Output | Send_Flag | Decision
%     ------------------|------------|-----------------|-----------|------------------
%            1 (Z)      | r <  Psend | Input pulse     |     1     | SEND
%            1 (Z)      | r >= Psend | Vacuum          |     0     | NOT_SEND
%            0 (X)      |    N/A     | Vacuum          |     0     | BLOCKED_X_WINDOW
%
% OptiSystem wiring:
%   Inputs  : Port 1 = Electrical (window selection: 1 = Z, 0 = X)
%             Port 2 = Optical    (coherent pulse train; mu already set
%                                  by upstream attenuation/modulation —
%                                  this component must not alter pulse
%                                  intensity on SEND)
%   Outputs : Port 1 = Optical    (SEND: passthrough, otherwise: vacuum)
%   Main tab: Sampled signal domain = Time
%   Run command: SNSDecisionBlock
%   User Parameters: Parameter0 = Psend (default 0.3, user-editable)
%                    Parameter1 = Mu (constant send intensity, default 0.1, user-editable)
%
% Modeling assumption:
%   Each element of InputPort1's electrical signal vector corresponds to
%   exactly one transmission window, and the optical signal on
%   InputPort2 carries exactly one sample per window per channel/
%   polarization component (a system-level abstraction consistent with
%   the Window Selection block this component receives its control
%   signal from). If your transmitter oversamples multiple points per
%   window, downsample/align before this block or tell me the
%   samples-per-window count and I'll adjust the indexing.
%
% Side effect:
%   Reads the existing alice_signal_record.csv (written by the Window
%   Selection block), merges in Send_Flag, Random_Number, Mu, and
%   Decision columns per Window_Index — preserving the existing
%   Selected_Window column — and rewrites the file.
%==========================================================================


outputDirectory = 'D:\drdo\osd\optisystem_prototype\matlab_cosim';
csvFileName     = fullfile(outputDirectory, 'alice_signal_record.csv');

if exist('Parameter0', 'var') && ~isempty(Parameter0)
    Psend = Parameter0;
else
    Psend = 0.3;
end

if exist('Parameter1', 'var') && ~isempty(Parameter1)
    Mu = Parameter1;
else
    Mu = 0.1;
end


if ~strcmp(InputPort1.TypeSignal, 'Electrical')
    error('SNSDecisionBlock:InvalidControlSignal', ...
        'InputPort1 must be the Electrical window-selection signal.');
end
if ~strcmp(InputPort2.TypeSignal, 'Optical')
    error('SNSDecisionBlock:InvalidOpticalSignal', ...
        'InputPort2 must be the Optical coherent-pulse signal.');
end

controlSignal = InputPort1.Sampled.Signal;   % 1xN: 1 = Z-window, 0 = X-window
numWindows    = length(controlSignal);


[ls, numChannels] = size(InputPort2.Sampled);
if ls == 0
    error('SNSDecisionBlock:EmptyOpticalSignal', ...
        'InputPort2.Sampled is empty — no optical channel to process.');
end


for c = 1:numChannels
    if size(InputPort2.Sampled(1, c).Signal, 2) ~= numWindows
        error('SNSDecisionBlock:SizeMismatch', ...
            ['Optical channel %d has %d samples but the electrical ' ...
             'control signal has %d windows. Align sample rates ' ...
             'before this block.'], c, ...
            size(InputPort2.Sampled(1, c).Signal, 2), numWindows);
    end
end


OutputPort1 = InputPort2;

sendFlagLog = zeros(numWindows, 1);
randomLog   = nan(numWindows, 1);
muLog       = zeros(numWindows, 1);
decisionLog = strings(numWindows, 1);

for w = 1:numWindows

    isZWindow = controlSignal(w) > 0.5;

    if isZWindow
        r = rand;

        if r < Psend
            sendFlagLog(w) = 1;
            randomLog(w)   = r;
            decisionLog(w) = "SEND";
            muLog(w)       = Mu;   % constant send intensity, not measured power

            for c = 1:numChannels
                pulseSample = InputPort2.Sampled(1, c).Signal(:, w);
                OutputPort1.Sampled(1, c).Signal(:, w) = pulseSample;
            end

        else
        
            sendFlagLog(w) = 0;
            randomLog(w)   = r;
            decisionLog(w) = "NOT_SEND";
            muLog(w)       = 0;  

            for c = 1:numChannels
                OutputPort1.Sampled(1, c).Signal(:, w) = 0;  % vacuum
            end
        end

    else
       
        sendFlagLog(w) = 0;
        randomLog(w)   = NaN;
        decisionLog(w) = "BLOCKED_X_WINDOW";
        muLog(w)       = 0;

        for c = 1:numChannels
            OutputPort1.Sampled(1, c).Signal(:, w) = 0;  % vacuum
        end
    end

end


windowIndex = (1:numWindows)';

newColumns = table(windowIndex, sendFlagLog, randomLog, muLog, decisionLog, ...
    'VariableNames', {'Window_Index', 'Send_Flag', 'Random_Number', 'Mu', 'Decision'});

if exist(csvFileName, 'file')
    existingTable = readtable(csvFileName);

    if height(existingTable) == numWindows
        mergedTable = join(existingTable, newColumns, 'Keys', 'Window_Index');
    else
    
        warning('SNSDecisionBlock:RowCountMismatch', ...
            ['Existing CSV has %d rows but this run processed %d ' ...
             'windows. Writing SNS results without merging ' ...
             'Selected_Window.'], height(existingTable), numWindows);
        mergedTable = newColumns;
    end
else
    warning('SNSDecisionBlock:MissingUpstreamLog', ...
        'alice_signal_record.csv not found at %s; creating it fresh.', ...
        csvFileName);
    mergedTable = newColumns;
end

writetable(mergedTable, csvFileName);