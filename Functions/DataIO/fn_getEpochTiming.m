function timing = fn_getEpochTiming(protocol, fs, motor_buffer_s, recovery_delay_s)
% FN_GETEPOCHTIMING  Build the Baseline/Evoked/Recovery epoch timing struct
% for a given Aplysia perturbation protocol.
%
%   timing = FN_GETEPOCHTIMING(protocol, fs, motor_buffer_s, recovery_delay_s)
%
% INPUTS
%   protocol         - '12min' or '20min'
%   fs               - sampling rate (fps)
%   motor_buffer_s   - buffer (s) added after P9 motor onset before calling
%                      the epoch "Evoked" (excludes motor transient)
%   recovery_delay_s - delay (s) added after C2 perturbation before calling
%                      the epoch "Recovery"
%
% OUTPUT
%   timing - struct with fields:
%     .t_P9, .t_C2, .t_end            - protocol landmark times (s)
%     .f_P9, .f_C2                    - landmark times in frames
%     .t_evoked_start, .t_recovery_start
%     .wins_s        - 3x2 [start end] (s) for Baseline/Evoked/Recovery
%     .win_labels    - {'Baseline','Evoked','Recovery'}

switch protocol
    case '12min'
        t_P9 = 120;  t_C2 = 420;  t_end = 1172880/fs;
    case '20min'
        t_P9 = 300;  t_C2 = 600;  t_end = 1953792/fs;
    otherwise
        error('fn_getEpochTiming:badProtocol', ...
            'Unknown protocol "%s" - expected ''12min'' or ''20min''.', protocol);
end

t_evoked_start   = t_P9 + motor_buffer_s;
t_recovery_start = t_C2 + recovery_delay_s;

timing.t_P9   = t_P9;
timing.t_C2   = t_C2;
timing.t_end  = t_end;
timing.f_P9   = round(t_P9 * fs);
timing.f_C2   = round(t_C2 * fs);
timing.t_evoked_start   = t_evoked_start;
timing.t_recovery_start = t_recovery_start;

timing.wins_s = [0,               t_P9; ...
                 t_evoked_start,  t_C2; ...
                 t_recovery_start, t_end];
timing.win_labels = {'Baseline', 'Evoked', 'Recovery'};

fprintf('Epoch definitions:\n');
fprintf('  Baseline : 0 - %.0f s\n', t_P9);
fprintf('  Evoked   : %.0f - %.0f s  (P9 + %ds buffer)\n', t_evoked_start, t_C2, motor_buffer_s);
fprintf('  Recovery : %.0f - %.0f s  (C2 + %ds delay)\n', t_recovery_start, t_end, recovery_delay_s);

end
