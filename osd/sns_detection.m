
OutputPort1 = InputPort1;
OutputPort2 = InputPort1;

epsilon = 0.5;  % Threshold for SEND decision (default: 0.5)
verbose = true;        % Enable debug output


if strcmp(InputPort1.TypeSignal, 'Optical')
   
    if isfield(InputPort2, 'Sampled') && ~isempty(InputPort2.Sampled)
        ctrl_samples = InputPort2.Sampled.Signal;
    else
        error('No electrical control signal on Input Port 2');
    end
    

    [rows, cols] = size(InputPort1.Sampled);
    
    if rows > 0
        for counter = 1:cols
            % Get optical samples for this signal
            opt_samples = InputPort1.Sampled(1, counter).Signal;
            
            % Get number of samples
            num_samples = size(opt_samples, 2);
            
            % Align control samples length if needed
            if length(ctrl_samples) ~= num_samples
                % Use the smaller length
                num_samples = min(length(ctrl_samples), num_samples);
                opt_samples = opt_samples(:, 1:num_samples);
                ctrl_trimmed = ctrl_samples(1:num_samples);
            else
                ctrl_trimmed = ctrl_samples;
            end
            
            % Handle complex control signal (use real part)
            if ~isreal(ctrl_trimmed)
                ctrl_values = real(ctrl_trimmed);
            else
                ctrl_values = ctrl_trimmed;
            end
            

            send_mask = ctrl_values >= epsilon;
            send_count = sum(send_mask);
            notsend_count = num_samples - send_count;
            
            if verbose
                fprintf('SNS Decision Switch:\n');
                fprintf('  Signal %d: %d samples\n', counter, num_samples);
                fprintf('  SEND: %d (%.1f%%)\n', send_count, 100*send_count/num_samples);
                fprintf('  NOT SEND: %d (%.1f%%)\n', notsend_count, 100*notsend_count/num_samples);
            end
          
            % SEND path: pass optical signal when mask is true
            send_output = zeros(size(opt_samples));
            send_output(:, send_mask) = opt_samples(:, send_mask);
            
            % NOT SEND path: vacuum (all zeros)
            notsend_output = zeros(size(opt_samples));
            
       
            OutputPort1.Sampled(1, counter).Signal = send_output;
            OutputPort2.Sampled(1, counter).Signal = notsend_output;
        end
    end
    

    if isfield(InputPort1, 'Parameterized') && ~isempty(InputPort1.Parameterized)
        % For parameterized signals, we need to handle power
        % Since this is a switch (not attenuator), we either pass or block
        % The decision is already applied above to sampled signals
        % Copy parameterized signals as-is (they don't carry bit-level info)
        OutputPort1.Parameterized = InputPort1.Parameterized;
        OutputPort2.Parameterized = InputPort1.Parameterized;
    end
    

    if isfield(InputPort1, 'Noise') && ~isempty(InputPort1.Noise)
        OutputPort1.Noise = InputPort1.Noise;
        OutputPort2.Noise = InputPort1.Noise;
    end
    
else
    error('InputPort1 is not an optical signal');
end