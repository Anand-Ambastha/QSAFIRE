%==========================================================================
% OptiSystem MATLAB Component: Binary Window Router
%--------------------------------------------------------------------------
% Purpose:
%   Reads a single electrical input signal (binary samples, one sample per
%   transmission window) and routes each window to two complementary
%   electrical outputs according to a fixed truth table. Every routing
%   decision is logged to a CSV file.
%
%   Truth table:
%       Input | Output1 | Output2 | Window label
%       ------|---------|---------|-------------
%         1   |    1    |    0    |     Z
%         0   |    0    |    1    |     X
%
% OptiSystem wiring:
%   Inputs  : 1 port, Signal type = Electrical
%   Outputs : 2 ports, Signal type = Electrical, Electrical
%   Main tab: Sampled signal domain = Time
%   Run command: AliceWindowRouter
%
% Variables OptiSystem expects when this script finishes:
%   OutputPort1, OutputPort2
%
% Side effect:
%   Overwrites 'alice_signal_record.csv' in the current working directory
%   at the start of every run, then writes one row per processed sample:
%   Window_Index, Selected_Window
%==========================================================================


outputDirectory = 'D:\drdo\osd\optisystem_prototype\matlab_cosim';

if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end

logFileName = fullfile(outputDirectory, 'alice_signal_record.csv');

decisionThreshold = 0.5;   % values above this are treated as binary '1'


if ~strcmp(InputPort1.TypeSignal, 'Electrical')
    error('AliceWindowRouter:InvalidSignalType', ...
        'This component expects an Electrical input signal on InputPort1.');
end

inputSignal = InputPort1.Sampled.Signal;   % 1xN vector, one sample per window
numSamples  = length(inputSignal);


output1      = zeros(1, numSamples);
output2      = zeros(1, numSamples);
windowLabels = strings(1, numSamples);     % holds "Z" or "X" per window


for idx = 1:numSamples

    bitIsOne = inputSignal(idx) > decisionThreshold;

    if bitIsOne
        % Case 1: Input = 1 -> Z-window
        output1(idx)      = 1;
        output2(idx)      = 0;
        windowLabels(idx) = "Z";
    else
        % Case 2: Input = 0 -> X-window
        output1(idx)      = 0;
        output2(idx)      = 1;
        windowLabels(idx) = "X";
    end

end


OutputPort1 = InputPort1;
OutputPort1.Sampled.Signal = output1;

OutputPort2 = InputPort1;
OutputPort2.Sampled.Signal = output2;

windowIndex    = (1:numSamples)';
selectedWindow = windowLabels(:);

logTable = table(windowIndex, selectedWindow, ...
    'VariableNames', {'Window_Index', 'Selected_Window'});


writetable(logTable, logFileName);