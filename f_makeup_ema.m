%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential moving average (EMA) Make-Up Gain algorithm
%
% Author: Jay Piamjariyakul
%
% Sources
% EMA difference equation:
% - D. Ward, "Applications of loudness models in audio engineering," Ph.D. 
%   dissertation, Birmingham City University, 2017
% K-weighting documentation:
% - International Telecommunications Union, "ITU-R BS.1770 Algorithms to 
%   measure audioprogramme loudness and true-peak audio level," 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [out_L_ema, timeTaken] = f_makeup_ema(in_refr, out_fil)

Fs_filter = 44100;
sz_channel = size(in_refr, 2);
%% Obtain RLB & Pre-K coefficients
[coef_rlb_b, coef_rlb_a] = f_getCoef_rlb(Fs_filter); 
[coef_prK_b, coef_prK_a] = f_getCoef_preK(Fs_filter); 

%% Set necessary parameters
% Time constant of EMA
const_ema = 0.125;
coef_ema = 1 - exp(-1 / (Fs_filter * const_ema));
% Initialises output memory
out_L_ema = zeros(size(in_refr, 1), sz_channel);

% Previous loudness values
mem_ema_L_i = zeros(1, sz_channel);
mem_ema_L_f = zeros(1, sz_channel);

% Memory of previous delay values in 1D filter
delay_f_prK = zeros(2, sz_channel);
delay_f_rlb = zeros(2, sz_channel);
delay_i_prK = zeros(2, sz_channel);
delay_i_rlb = zeros(2, sz_channel);
% Stores buffer for filters

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

%%% 3) Finds EMA sequence & stores them
% Finds EMA of filtr
mem_ema_L_f = (coef_ema * (val_proc_filtr .^ 2)) + ((1 - coef_ema) * mem_ema_L_f);
% Finds EMA of input
mem_ema_L_i = (coef_ema * (val_proc_input .^ 2)) + ((1 - coef_ema) * mem_ema_L_i);
% Converts vals to log form
log_L_filtr = 10 * log10(mem_ema_L_f);
log_L_input = 10 * log10(mem_ema_L_i);

% Standard "ITU-R BS.1770-3" requires subtracting these LOG vals with –0.691
% However, we intend to subtract the values together, so this becomes redundant

%%% 4) Calculate the correcting gain
% Find the error value
err_L = log_L_input - log_L_filtr;
% Converts err value to gain
gain_corrector = 10 .^ (err_L ./ 20);

%%% 5) Corrects value
% Corrects NaN values - set to zero
gain_corrector(isnan(gain_corrector)) = 0;
% Prevents infinite gain - set to max floating point
gain_corrector(isinf(gain_corrector)) = sign(gain_corrector(isinf(gain_corrector))) * realmax;

%%% 6) Apply gain to delayed filter output at sample n
out_L_ema(i_sample, :) = out_fil(i_sample, :) .* gain_corrector; % Get current output
    
end
timeTaken = toc;

end