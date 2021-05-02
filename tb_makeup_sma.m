%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple moving average (SMA) Make-Up Gain algorithm
%
% Author: Jay Piamjariyakul
%
% Sources
% RMS summation (SMA was derived from the summation):
% - G. A. Soulodre, "Evaluation of objective loudness meters," in Audio 
%   Engineering SocietyConvention 116, Audio Engineering Society, 2004.
% K-weighting documentation:
% - International Telecommunications Union, "ITU-R BS.1770 Algorithms to 
%   measure audioprogramme loudness and true-peak audio level," 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [out_L_sma, timeTaken] = tb_makeup_sma(in_refr, out_fil)

Fs_filter = 44100;
sz_channel = size(in_refr, 2);
%% Obtain RLB & Pre-K coefficients
[coef_rlb_b, coef_rlb_a] = f_getCoef_rlb(Fs_filter); 
[coef_prK_b, coef_prK_a] = f_getCoef_preK(Fs_filter); 

%% Set necessary parameters
size_window = round(Fs_filter * 0.4); % 0.4s from ITU document
% size_window = round(Fs_filter * 0.4); % 0.4s from ITU document
% Initialises output memory
out_L_sma = zeros(size(in_refr, 1), sz_channel);

% Memory of x[n] and x[n-W]
mem_K_f = zeros(size_window, sz_channel);
mem_K_i = zeros(size_window, sz_channel);

% Previous loudness values
mem_sma_L_i = zeros(1, sz_channel);
mem_sma_L_f = zeros(1, sz_channel);

% Memory of previous delay values in 1D filter
delay_f_prK = zeros(2, sz_channel); % Stores buffer for filtr algorithm
delay_f_rlb = zeros(2, sz_channel); % Stores buffer for filtr algorithm
delay_i_prK = zeros(2, sz_channel); % Stores buffer for filtr algorithm
delay_i_rlb = zeros(2, sz_channel); % Stores buffer for filtr algorithm

%% Run make-up gain
tic;
for i_sample = 1:length(in_refr)

%%% 1) Obtain sequence of signals to process
val_proc_filtr = out_fil(i_sample, :);
val_proc_input = in_refr(i_sample, :);

%%% 2) Filter signals
% Applies K weighting on filtered signal
[val_proc_filtr, delay_f_prK] = f_1dFilter(coef_prK_b, coef_prK_a, val_proc_filtr, delay_f_prK);
[val_proc_filtr, delay_f_rlb] = f_1dFilter(coef_rlb_b, coef_rlb_a, val_proc_filtr, delay_f_rlb);
% Applies K weighting on reference signal
[val_proc_input, delay_i_prK] = f_1dFilter(coef_prK_b, coef_prK_a, val_proc_input, delay_i_prK);
[val_proc_input, delay_i_rlb] = f_1dFilter(coef_rlb_b, coef_rlb_a, val_proc_input, delay_i_rlb);
% Processes the first and final values of 'window'
mem_K_f = [mem_K_f(2:end, :); val_proc_filtr];
mem_K_i = [mem_K_i(2:end, :); val_proc_input];

%%% 3) Calculate Leq values
mem_sma_L_f = mem_sma_L_f + (mem_K_f(end, :).^2 - mem_K_f(1, :).^2)/size_window;
mem_sma_L_i = mem_sma_L_i + (mem_K_i(end, :).^2 - mem_K_i(1, :).^2)/size_window;
Leq_filtr = mem_sma_L_f;
Leq_input = mem_sma_L_i;
% Converts linear value to LOG
log_L_filtr = 10 * log10(Leq_filtr);
log_L_input = 10 * log10(Leq_input);

% Standard "ITU-R BS.1770-3" requires subtracting these LOG vals with –0.691
% However, we intend to subtract the values together, so this becomes redundant

%%% 4) Calculate the correcting gain
% Find the error value
err_L = log_L_input - log_L_filtr;
% Converts err value to gain
gain_corrector = real(10 .^ (err_L ./ 20));

%%% 5) Corrects value
% Corrects NaN values - set to zero
gain_corrector(isnan(gain_corrector)) = 0;
% Prevents infinite gain - set to max floating point
gain_corrector(isinf(gain_corrector)) = sign(gain_corrector(isinf(gain_corrector))) * realmax;

%%% 6) Apply gain to delayed filter output at sample n
out_now = out_fil(i_sample, :) .* gain_corrector; % Get current output
out_L_sma(i_sample, :) = out_now;

end
timeTaken = toc;

end