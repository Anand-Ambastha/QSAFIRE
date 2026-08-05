
%==========================================================================
% Charlie X-Basis Detection Event Logger
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
% Output Files:
%   D:\drdo\osd\code\csv\charlie_record_x.csv
%   D:\drdo\osd\code\csv\charlie_xbasis_log.csv
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

FullLogFileName    = fullfile(outputDirectory,'charlie_record_x.csv');
SummaryLogFileName = fullfile(outputDirectory,'charlie_xbasis_log.csv');
WindowSelectFile   = fullfile(outputDirectory,'window_selection_log.csv');


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


if exist(WindowSelectFile,'file')

    fid_chk = fopen(WindowSelectFile,'r');

    if fid_chk ~= -1

        fgetl(fid_chk);

        lineCount = 0;

        while ~feof(fid_chk)

            line = fgetl(fid_chk);

            if ischar(line) && ~isempty(strtrim(line))
                lineCount = lineCount + 1;
            end

        end

        fclose(fid_chk);

        if lineCount ~= NumWindows

            warning(['Window count mismatch: ' ...
                     'window_selection_log.csv = %d, ' ...
                     'Charlie X windows = %d'], ...
                     lineCount, NumWindows);

        end

    end

end


ClickFlag = zeros(1,NumWindows);

-
fidFull = fopen(FullLogFileName,'w');

if fidFull == -1
    error('Cannot open %s',FullLogFileName);
end

fprintf(fidFull,...
'Window_Index,SPD1,SPD2,Detector,Event\n');

fidSummary = fopen(SummaryLogFileName,'w');

if fidSummary == -1

    fclose(fidFull);
    error('Cannot open %s',SummaryLogFileName);

end

fprintf(fidSummary,...
'Window_Index,Detector,Event\n');

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

    fprintf(fidFull,...
        '%d,%d,%d,%s,%s\n',...
        WindowIndex,...
        SPD1_bit,...
        SPD2_bit,...
        Detector,...
        Event);

    fprintf(fidSummary,...
        '%d,%s,%s\n',...
        WindowIndex,...
        Detector,...
        Event);

    ClickFlag(w) = double(click1 || click2);

end


fclose(fidFull);
fclose(fidSummary);


OutputPort1.Sampled.Signal = ClickFlag;

if isfield(InputPort1.Sampled,'Time')
    OutputPort1.Sampled.Time = InputPort1.Sampled.Time;
end


disp(['Charlie X Logger: Processed ',num2str(NumWindows),' windows']);
disp(['Saved: ',FullLogFileName]);
disp(['Saved: ',SummaryLogFileName]);
