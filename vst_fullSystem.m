%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make-Up Gain and Moog VCF audio plugin
%
% Author: Jay Piamjariyakul
%
% Sources
% Nonlinear digital Moog VCF model:
% - P. Daly, "A comparison of virtual analogue Moog VCF models," Master's
%   thesis, Univ. ofEdinburgh, Edinburgh, UK, Aug, 2012.
% EMA difference equation:
% - D. Ward, "Applications of loudness models in audio engineering," Ph.D. 
%   dissertation, Birmingham City University, 2017
% K-weighting documentation:
% - International Telecommunications Union, "ITU-R BS.1770 Algorithms to 
%   measure audioprogramme loudness and true-peak audio level," 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef vst_fullSystem < audioPlugin % Inherits audioPlugin
properties % Tunable properties
    k = 3.99; % Feedback gain
    Fc = 440;   % Cutoff frequency
    factor_gain = 0; % Gain-scaling factor
end
properties (Access=private)        
    %% Variables used in K-weighting
    % Memory of previous delay values in 1D k-weighting filter
    delay_f_prK = zeros(2, 2);
    delay_f_rlb = zeros(2, 2);
    delay_i_prK = zeros(2, 2);
    delay_i_rlb = zeros(2, 2);
    % Stores buffer for filtr algorithm

    % Coefficients used in K-weighting
    coef_rlb_b = zeros(1, 3);
    coef_rlb_a = zeros(1, 3);
    coef_prK_b = zeros(1, 3);
    coef_prK_a = zeros(1, 3);

    %% Variables used in the VCF
    y_vcf = zeros(4, 2); % Stores the outputs of each stage
    % Ts = 1/Fs_filter; % 1/Fs == Ts: recall n=t/Ts == t*Fs
    y_vcf_d2 = zeros(1, 2); % Twice-delayed 

    %% Variables used in EMA loudness calculations
    % Previous loudness values
    mem_ema_L_i = zeros(1, 2);
    mem_ema_L_f = zeros(1, 2);        
    % EMA weighting coefficient
    coef_ema = 0;
end
    
properties (Constant)
%% Sets sampling frequency of system
Fs_filter = 44100;
const_ema = 0.125;
%sz_channel = 2;

% Map tunable property to plugin parameter.
PluginInterface = audioPluginInterface( ...
    audioPluginParameter(   'k', ...
                            'Mapping', {'lin', 0, 4}, ...
                            'Style', 'rotary', ...
                            'Layout', [2, 1], ...
                            'DisplayName', 'Feedback', ...
                            'DisplayNameLocation', 'Above' ...
                            ), ...
    audioPluginParameter(   'Fc', ...
                            'Label','Hz', ...
                            'Mapping', {'log', 1, 22000}, ...
                            'Style', 'rotary', ...
                            'Layout', [2, 2], ...
                            'DisplayName', 'Cutoff Frequency', ...
                            'DisplayNameLocation', 'Above' ...
                            ), ...
    audioPluginParameter(   'factor_gain', ...
                            'Mapping', {'lin', 0, 100}, ...
                            'Style', 'rotary', ...
                            'Layout', [4, 2], ...
                            'Label', '%', ...
                            'DisplayName', 'Make-Up Gain Strength', ...
                            'DisplayNameLocation', 'Above' ...
                            ), ...
    audioPluginGridLayout( ... 
        'RowHeight',[20, 100, 20, 100], ... 
        'ColumnWidth',[150, 150] ...
        ) ...
    );

end

methods
function out_L_rtp = process(plugin, in_refr) % Audio processing methods
    %% Other variables (dynamic!)
    g = f_coefWeight(plugin.Fc, plugin.Fs_filter);

    %% Set necessary parameters
    % Initialises output memory
    out_L_rtp = zeros(size(in_refr, 1), 2);

    %% Run VCF and make-up gain
    for i_sample = 1:length(in_refr)
        %%% 0) VCF module: Filter current value
        y_vcf_0 = in_refr(i_sample, :) - ( plugin.k * ( 0.5 * ( plugin.y_vcf(4, :) + plugin.y_vcf_d2 ) ) );
        plugin.y_vcf(1, :) = plugin.y_vcf(1, :) + ( g * ( tanh( y_vcf_0     ) - tanh( plugin.y_vcf(1, :) ) ) );
        plugin.y_vcf(2, :) = plugin.y_vcf(2, :) + ( g * ( tanh( plugin.y_vcf(1, :) ) - tanh( plugin.y_vcf(2, :) ) ) );
        plugin.y_vcf(3, :) = plugin.y_vcf(3, :) + ( g * ( tanh( plugin.y_vcf(2, :) ) - tanh( plugin.y_vcf(3, :) ) ) );
        plugin.y_vcf_d2 = plugin.y_vcf(4, :); % Stores twice-delay output before overriding
        plugin.y_vcf(4, :) = plugin.y_vcf(4, :) + ( g * ( tanh( plugin.y_vcf(3, :) ) - tanh( plugin.y_vcf(4, :) ) ) );

        %%% 1) Obtain sequence of signals to process
        out_L_rtp(i_sample, :) = plugin.y_vcf(4, :);
        val_proc_filtr = out_L_rtp(i_sample, :);
        val_proc_input = in_refr(i_sample, :);

        %%% 2) Filter signals
        % Applies K weighting on filtered signal
        [val_proc_filtr, plugin.delay_f_prK] = f_filter1D(plugin.coef_prK_b, plugin.coef_prK_a, val_proc_filtr, plugin.delay_f_prK);
        [val_proc_filtr, plugin.delay_f_rlb] = f_filter1D(plugin.coef_rlb_b, plugin.coef_rlb_a, val_proc_filtr, plugin.delay_f_rlb);
        % Applies K weighting on reference signal
        [val_proc_input, plugin.delay_i_prK] = f_filter1D(plugin.coef_prK_b, plugin.coef_prK_a, val_proc_input, plugin.delay_i_prK);
        [val_proc_input, plugin.delay_i_rlb] = f_filter1D(plugin.coef_rlb_b, plugin.coef_rlb_a, val_proc_input, plugin.delay_i_rlb);

        %%% 3) Finds EMA sequence & stores them
        % Finds EMA of filtr
        plugin.mem_ema_L_f = (plugin.coef_ema * (val_proc_filtr .^ 2)) + ((1 - plugin.coef_ema) * plugin.mem_ema_L_f);
        % Finds EMA of input
        plugin.mem_ema_L_i = (plugin.coef_ema * (val_proc_input .^ 2)) + ((1 - plugin.coef_ema) * plugin.mem_ema_L_i);
        % Converts vals to log form
        log_L_filtr = log10(plugin.mem_ema_L_f);
        log_L_input = log10(plugin.mem_ema_L_i);

        % Standard "ITU-R BS.1770-3" requires subtracting these LOG vals with –0.691
        % However, we intend to subtract the values together, so this becomes redundant

        %%% 4) Calculate the correcting gain
        % Find the error value
        err_L = log_L_input - log_L_filtr;
        % Converts err value to gain
        gain_corrector = 10 .^ (err_L ./ 2);

        %%% 5) Corrects value
        % Corrects NaN values - set to zero
        gain_corrector(isnan(gain_corrector)) = 0;

        %%% 6) Apply gain to delayed filter output at sample n
        out_L_rtp(i_sample, :) = 0.01 * (((100 - plugin.factor_gain) .* plugin.y_vcf(4, :)) + (plugin.factor_gain .* (plugin.y_vcf(4, :) .* gain_corrector))); % Get current output          
        % The gain strength controls how much impact does the gain stage have
        % IE if the gain strength is 0% then the system outputs only the VCF output
        % If the gain strength is 100% then the system outputs the gain-applied output wholly
        % Gain strength at 50% outputs half the strength of the make-up gain.

    end % END AGC      

    %% Applies 1D filter
    function [out_filt, seq_delay] = f_filter1D(b, a, in_audio, seq_delay)
        b0 = b(1); b1 = b(2); b2 = b(3);
        a1 = a(2); a2 = a(3);

        %plugin.sz_channel = size(in_audio, 2);
        out_filt = zeros(size(in_audio, 1), 2);
        for i_filter = 1:size(in_audio, 1)
            delayNow = in_audio(i_filter, :) - (a1 * seq_delay(1, :)) - (a2 * seq_delay(2, :));
            out_filt(i_filter, :) = (b0 * delayNow) + (b1 * seq_delay(1, :)) + (b2 * seq_delay(2, :));
            seq_delay(2, :) = seq_delay(1, :);
            seq_delay(1, :) = delayNow;
        end
    end
    %% Obtains weighting coefficient
    function coefWeight = f_coefWeight(Fc, Fs)
        coefWeight = 1 - exp( (-2 * pi) * ( Fc / Fs ) );
    end
end

%% reset() function - runs at startup or manual reset
function reset(plugin)
    %%% Define channel size and sampling frequency
    plugin.delay_f_prK = zeros(2, 2);
    plugin.delay_f_rlb = zeros(2, 2);
    plugin.delay_i_prK = zeros(2, 2);
    plugin.delay_i_rlb = zeros(2, 2);
    % Stores buffer for filtr algorithm

    % Previous loudness values
    plugin.mem_ema_L_i = zeros(1, 2);
    plugin.mem_ema_L_f  = zeros(1, 2);

    %% Variables used in the VCF
    plugin.y_vcf = zeros(4, 2);
    % Ts = 1/Fs_filter; % 1/Fs == Ts: recall n=t/Ts == t*Fs
    plugin.y_vcf_d2 = zeros(1, 2); % Twice-delayed 

    %% Obtain RLB & Pre-K coefficients
    [plugin.coef_rlb_b, plugin.coef_rlb_a] = f_getCoef_rlb(plugin.Fs_filter); 
    [plugin.coef_prK_b, plugin.coef_prK_a] = f_getCoef_preK(plugin.Fs_filter); 

    % Obtains the EMA coefficient
    plugin.coef_ema = 1 - exp(-1 / (plugin.Fs_filter * plugin.const_ema));

    %% Obtains pre-K filter's coefficients
    function [fil_preK_coef_b, fil_preK_coef_a] = f_getCoef_preK(Fs_filter)
        fil_preK_fc = 1681.9744510;
        fil_preK_q = 0.7071752;
        fil_preK_VL = 1;
        fil_preK_VB = 1.2587209;
        fil_preK_VH = 1.5848647;
        fil_preK_Omega = tan(pi * (fil_preK_fc / Fs_filter)); % The only fs-dependent parameter

        fil_preK_coef_a = [  (fil_preK_Omega ^ 2) + (fil_preK_Omega / fil_preK_q) + 1,  ... % Coefficient a0
                            2 * ((fil_preK_Omega^2) - 1),    ...                            % Coefficient a1
                            (fil_preK_Omega^2) - (fil_preK_Omega / fil_preK_q) + 1 ...      % Coefficient a2
                            ];
        fil_preK_coef_b = [  (fil_preK_VL * fil_preK_Omega^2) + (fil_preK_VB * fil_preK_Omega / fil_preK_q) + fil_preK_VH, ...  % Coefficient b0
                            2 * ((fil_preK_VL * fil_preK_Omega^2) - fil_preK_VH),  ...                                          % Coefficient b1
                            (fil_preK_VL * fil_preK_Omega^2) - (fil_preK_VB * fil_preK_Omega / fil_preK_q) + fil_preK_VH  ...   % Coefficient b2
                            ];
        fil_preK_coef_a0 = fil_preK_coef_a(1);
        fil_preK_coef_a = fil_preK_coef_a / fil_preK_coef_a0; % Normalise coefficients to a0
        fil_preK_coef_b = fil_preK_coef_b / fil_preK_coef_a0; % Normalise coefficients to a0
    end

    %% Obtains RLB filter's coefficients
    function [fil_rlb_coef_b, fil_rlb_coef_a] = f_getCoef_rlb(Fs_filter)
        fil_rlb_fc = 38.1354709;
        fil_rlb_q = 0.5003270;
        fil_rlb_VL = 0;
        fil_rlb_VB = 0;
        fil_rlb_VH = 1.0049949;
        fil_rlb_Omega = tan(pi * (fil_rlb_fc / Fs_filter)); % The only fs-dependent param

        fil_rlb_coef_a = [  (fil_rlb_Omega ^ 2) + (fil_rlb_Omega / fil_rlb_q) + 1,  ... % Coefficient a0
                            2 * ((fil_rlb_Omega^2) - 1),    ...                         % Coefficient a1
                            (fil_rlb_Omega^2) - (fil_rlb_Omega / fil_rlb_q) + 1 ...     % Coefficient a2
                            ];
        fil_rlb_coef_b = [  (fil_rlb_VL * fil_rlb_Omega^2) + (fil_rlb_VB * fil_rlb_Omega / fil_rlb_q) + fil_rlb_VH, ... % Coefficient b0
                            2 * ((fil_rlb_VL * fil_rlb_Omega^2) - fil_rlb_VH),  ...                                     % Coefficient b1
                            (fil_rlb_VL * fil_rlb_Omega^2) - (fil_rlb_VB * fil_rlb_Omega / fil_rlb_q) + fil_rlb_VH  ... % Coefficient b2
                            ];
        fil_rlb_coef_a0 = fil_rlb_coef_a(1);
        fil_rlb_coef_b = fil_rlb_coef_b / fil_rlb_coef_a0; % Normalise coefficients to a0
        fil_rlb_coef_a = fil_rlb_coef_a / fil_rlb_coef_a0; % Normalise coefficients to a0
    end

end      
    end % END PROCESS
end