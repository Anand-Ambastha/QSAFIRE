
%==========================================================================
% Charlie Z-Basis Detection Event Logger
% SNS-TF-QKD
%
% One sample = one transmission window
%
% Inputs:
%   InputPort1 = SPD1 electrical output
%   InputPort2 = SPD2 electrical output
%
% Outputs:
%   OutputPort1 = Click flag (1 if any detector clicked, else 0)
%
% User Parameters:
%   Parameter0 = Threshold
%   Parameter1 = StartWindowIndex
%
% Output File:
%   D:\drdo\osd\code\csv\charlie_record_z.csv
%
% Detection Logic:
%
% SPD1 SPD2 Result
% ----------------------
%  1    0    D1_CLICK
%  0    1    D2_CLICK
%  1    1    DOUBLE_CLICK
%  0    0    NO_CLICK
%
%==========================================================================

OutputPort1 = InputPort1;


if exist('Parameter0','var') && ~isempty(Parameter0)
    Threshold = Parameter0;
else
    Threshold = 1.2e-9;
end

if exist('Parameter1','var') && ~isempty(Parameter1)
    StartWindowIndex = Parameter1;
else
    StartWindowIndex = 1;
end

outputDirectory = 'D:\drdo\osd\optisystem_prototype\outputs';

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

OutputFileName = fullfile(outputDirectory,'charlie_record_z.csv');

if ~strcmp(InputPort1.TypeSignal,'Electrical')
    error('InputPort1 must be an Electrical signal.');
end

if ~strcmp(InputPort2.TypeSignal,'Electrical')
    error('InputPort2 must be an Electrical signal.');
end


Signal1 = InputPort1.Sampled.Signal;
Signal2 = InputPort2.Sampled.Signal;

if length(Signal1) ~= length(Signal2)
    error('SPD1 and SPD2 lengths do not match.');
end


NumWindows = length(Signal1);


ClickFlag = zeros(1,NumWindows);


fid = fopen(OutputFileName,'w');

if fid == -1
    error('Cannot open %s',OutputFileName);
end

fprintf(fid,...
'Window_Index,SPD1,SPD2,Detector,Event\n');

for w = 1:NumWindows

    click1 = Signal1(w) > Threshold;
    click2 = Signal2(w) > Threshold;

    SPD1_bit = double(click1);
    SPD2_bit = double(click2);


    if click1 && ~click2

        Detector = 'D1';
        Event    = 'D1_CLICK';

    elseif click2 && ~click1

        Detector = 'D2';
        Event    = 'D2_CLICK';

    elseif click1 && click2

        Detector = 'BOTH';
        Event    = 'DOUBLE_CLICK';

    else

        Detector = 'NONE';
        Event    = 'NO_CLICK';

    end

    WindowIndex = StartWindowIndex + w - 1;


    fprintf(fid,...
        '%d,%d,%d,%s,%s\n',...
        WindowIndex,...
        SPD1_bit,...
        SPD2_bit,...
        Detector,...
        Event);

    ClickFlag(w) = double(click1 || click2);

end

fclose(fid);

OutputPort1.Sampled.Signal = ClickFlag;

if isfield(InputPort1.Sampled,'Time')
    OutputPort1.Sampled.Time = InputPort1.Sampled.Time;
end


disp(['Charlie Z Logger: Processed ',num2str(NumWindows),' windows']);
disp(['Saved: ',OutputFileName]);

