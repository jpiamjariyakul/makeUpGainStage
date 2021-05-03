%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testbench
%
% Author: Jay Piamjariyakul
%
% Sources
% Dimensionless nonlinear digital Moog VCF difference equations:
% - P. Daly, "A comparison of virtual analogue Moog VCF models," Master's
%   thesis, Univ. ofEdinburgh, Edinburgh, UK, Aug, 2012.
% Original nonlinear VCF implementation:
% - A. Huovilainen, "Non-linear digital implementation of the Moog ladder
%   filter," in Proceed-ings of the International Conference on Digital 
%   Audio Effects (DAFx-04), 2004
% Original implementation of the RLB filter:
% - G. A. Soulodre, "Evaluation of objective loudness meters," in Audio 
%   Engineering SocietyConvention 116, Audio Engineering Society, 2004.
% Method for obtaining parameter values at any fs:
% - D. Ward, "Applications of loudness models in audio engineering," Ph.D. 
%   dissertation, Birmingham City University, 2017
% Original parameter values for fs=48000 Hz:
% - International Telecommunications Union, "ITU-R BS.1770 Algorithms to 
%   measure audioprogramme loudness and true-peak audio level," 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
clc, clear

%% 0) Defines essential parameters
Fs = 44100; % Sampling frequency (Hz)
Fc = 440;   % Cutoff frequency of Moog filter (Hz)
gainFeedbk = 3.99; % Feedback gain of Moog filter

range_f = (1 : Fs/2); % Frequency range of plots (half since mirror) (Hz)
range_f_norm = range_f * (2 * pi) / Fs;
duration = 1; % Signal duration (seconds)
t = 0:1/Fs:duration; t = t(2:end); % Time array (seconds)

% For plotting - define colors
dict_color = containers.Map({'blue', 'orange', 'yellow', 'purple', 'green', 'aqua', 'red'}, ...
                            {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'} ...
                            );

%% 1) Obtains test sequences
% Change 'flag' to switch between different test data
flag = "square";
switch flag
    case "impulse"
        in_audio = zeros(length(t), 2); in_audio(1, :) = 1;
    case "square"
        frqc_square = 220; % Change 'frqc_square' to change sq-wave frequency
        in_audio = square(2*pi * frqc_square * t); in_audio = [in_audio', in_audio'];
    case "audioFile" % Change 'filename' to desired audio file!
        filename = "./audio_testdata/waterNight_full.flac";
        [in_audio, Fs_now] = audioread(filename);
        if Fs_now ~= Fs % Resamples signal to sampling frequency specified
            in_audio = resample(in_audio, Fs, Fs_now);
        end
        t = 0:1/Fs:size(in_audio, 1)/Fs; t = t(2:end);
end

%% 2) Runs filter on specified signal
%Fc = 440; % Cutoff frequency (Hz)
% coefWeight_vcf = f_vcf_coefWeight(440, 44100);
coefWeight_vcf = 1 - exp( (-2 * pi) * ( Fc / Fs ) );
% gainFeedbk = 3.99;
out_filtr = f_runVcf(in_audio, coefWeight_vcf, gainFeedbk);
% Feedback gain = 3.99
% VCF weighting coefficient for cutoff frequency = 440 Hz

%% 3) Plots time response
% Specifies how long to plot the duration for
factor_t = length(t) * size(in_audio, 1) / length(t);
% Specifies snippet for subplot
t_subplt = zeros(1, 2);
t_subplt(1) = round(Fs * 0.055); % Plots from 55ms...
t_subplt(2) = round(Fs * 0.065); % ... to 65ms

figure;
hold on;
plot(t(1:factor_t), out_filtr(1:factor_t, 1));
title("Impulse Response of VCF Output (f_c = 440 kHz)");
xlabel("Duration (seconds)");
ylabel("Sample Amplitude (arb. units)");
grid on;
axes('position', [0.55 0.575 0.325 0.25]); box on;
hold on;
plot(t(t_subplt(1):t_subplt(2)), out_filtr(t_subplt(1):t_subplt(2)));
axis tight; grid on;

%% 4) Runs make-up gain on filtered signal
[out_L_sma, ~] = f_makeup_sma(in_audio, out_filtr);
[out_L_ema, ~] = f_makeup_ema(in_audio, out_filtr);

%% 5) Plots gain-applied outout
figure;
hold on;
    plot(t(1:factor_t), in_audio(1:factor_t, 1), "Color", dict_color("green"));
    plot(t(1:factor_t), out_L_sma(1:factor_t, 1), "Color", dict_color("green"));
    plot(t(1:factor_t), out_L_ema(1:factor_t, 1), "Color", dict_color("orange"));
    plot(t(1:factor_t), out_filtr(1:factor_t, 1), "Color", dict_color("blue"));
hold off;
title("Output of Gain Stage (VCF Impulse Response, f_c = 440 Hz)");
xlabel("Duration (seconds)");
ylabel("Sample Amplitude (arb. units)");
grid on;
legend(["SMA (RMS)", "EMA", "No make-up gain"], 'Location', 'southeast');

% Plots subplot
axes('position', [0.55 0.575 0.325 0.25]); box on; grid on;
hold on;
    plot(t(t_subplt(1):t_subplt(2)), in_audio(t_subplt(1):t_subplt(2)), "Color", dict_color("green"));
    plot(t(t_subplt(1):t_subplt(2)), out_L_sma(t_subplt(1):t_subplt(2)), "Color", dict_color("green"));
    plot(t(t_subplt(1):t_subplt(2)), out_L_ema(t_subplt(1):t_subplt(2)), "Color", dict_color("orange"));
    plot(t(t_subplt(1):t_subplt(2)), out_filtr(t_subplt(1):t_subplt(2)), "Color", dict_color("blue"));
hold off;
axis tight;

%% Plots gain-applied outout and other signals
% Plots input (reference) signal
figure;
subplot(3, 1, 1);
plot(t(1:factor_t), in_audio(1:factor_t, 1), "Color", dict_color("yellow"), "LineWidth", 1);
% ylim([-1.5, +1.5])
grid on; grid(gca,'minor');
subtitle("Original reference input (f = 220 Hz)", 'FontSize', 9);
title("Square Wave: Pre/Post Filter and Make-Up Gain Stages")

% Plots VCF output signal
subplot(3, 1, 2);
plot(t(1:factor_t), out_filtr(1:factor_t, 1), "Color", dict_color("blue"));
ylabel("Signal Amplitude (arb. units)");
% ylim([-0.5, +0.5])
grid on; grid(gca,'minor');
subtitle("VCF filtered output (f_c = 440 Hz)", 'FontSize', 9);

% Plots gain-applied output signals
subplot(3, 1, 3);
hold on;
plot(t(1:factor_t), out_L_sma(1:factor_t, 1), "Color", dict_color("green"));
plot(t(1:factor_t), out_L_ema(1:factor_t, 1), "Color", dict_color("orange"));
hold off;
xlabel("Duration (seconds)");
% ylim([-2, +2.5])
grid on; grid(gca,'minor');
subtitle("Make-up gain output", 'FontSize', 9);

%% Obtains frequency responses
% Obtains magnitude and phase responses
[m_in_audio, p_in_audio] = deal(abs(fft(in_audio)), angle(fft(in_audio)));
[m_out_filtr, p_out_filtr] = deal(abs(fft(out_filtr)), angle(fft(out_filtr)));
[m_out_L_sma, p_out_L_sma] = deal(abs(fft(out_L_sma)), angle(fft(out_L_sma)));
[m_out_L_ema, p_out_L_ema] = deal(abs(fft(out_L_ema)), angle(fft(out_L_ema)));

% Converts values to more useful forms
m_in_audio  = 10 * log( m_in_audio(range_f) );
p_in_audio  = mod(( unwrap(p_in_audio(range_f), 2 * pi) * 180/pi ), -360);
m_out_filtr = 10 * log( m_out_filtr(range_f) );
p_out_filtr = mod(( unwrap(p_out_filtr(range_f), 2 * pi) * 180/pi ), -360);
m_out_L_ema = 10 * log( m_out_L_ema(range_f) );
p_out_L_ema = mod(( unwrap(p_out_L_ema(range_f), 2 * pi) * 180/pi ), -360);
m_out_L_sma = 10 * log( m_out_L_sma(range_f) );
p_out_L_sma = mod(( unwrap(p_out_L_sma(range_f), 2 * pi) * 180/pi ), -360);

%% Plot frequency responses
figure;

% Plots magnitude response
subplot(2, 1, 1);
xlim([0, Fs/2]); ylim([-60, inf]);
hold on;
    plot(range_f, m_in_audio, "LineWidth", 0.5, "Color", dict_color("yellow"));
    plot(range_f, m_out_filtr, "LineWidth", 0.5, "Color", dict_color("blue"));
    plot(range_f, m_out_L_sma, "LineWidth", 1.0, "Color", dict_color("green"));
    plot(range_f, m_out_L_ema, "LineWidth", 1.0, "Color", dict_color("orange"));
hold off;
grid on; set(gca, 'XScale', 'log')
title("Magnitude Response of Make-Up Gain Stage (f_c = 440 Hz)");
ylabel("Magnitude (dB)"); xlabel("Frequency (Hz)");

% Plots phase response
subplot(2, 1, 2);
xlim([0, Fs/2]); ylim([-400, 45]);
hold on;
    plot(range_f, p_in_audio, "LineWidth", 0.5, "Color", dict_color("yellow"));
    plot(range_f, p_out_filtr, "LineWidth", 0.5, "Color", dict_color("blue"));
    plot(range_f, p_out_L_sma, "LineWidth", 1.0, "Color", dict_color("green"));
    plot(range_f, p_out_L_ema, "LineWidth", 1.0, "Color", dict_color("orange"));
hold off;
title("Phase Response of Make-Up Gain Stage (f_c = 440 Hz)");
ylabel("Phase /degrees"); xlabel("Frequency (Hz)");
grid on; set(gca, 'XScale', 'log'); set(gca, 'YTick', (-360:90:90));

legend(["No make-up gain", "SMA (RMS)", "EMA"], 'Location', 'southwest');

%% Compare the loudness between original, filtered, and gain-applied
sz_channel = size(in_audio, 2);

% Applies Ward's moving avg filtr
L_in_origin = zeros(length(in_audio), sz_channel);
L_out_filtr = zeros(length(in_audio), sz_channel);
L_out_L_sma = zeros(length(in_audio), sz_channel);
L_out_L_ema = zeros(length(in_audio), sz_channel);

% Defines EMA constant
const_loud = 0.125;
coef_loud = 1 - exp(-1 / (Fs * const_loud));
% Instantiates memory
mem_L_in_origin = zeros(1, sz_channel);
mem_L_out_filtr = zeros(1, sz_channel);
mem_L_out_L_sma = zeros(1, sz_channel);
mem_L_out_L_ema = zeros(1, sz_channel);

%%% EMA: Finds EMA sequence & stores them
for i_sample = 1:length(in_audio)
    % Get loudness of input reference
    mem_L_in_origin = (coef_loud * (in_audio(i_sample, :) .^ 2)) + ((1 - coef_loud) * mem_L_in_origin);
    L_in_origin(i_sample, :) = mem_L_in_origin;
    
    % Get loudness of VCF output
    mem_L_out_filtr = (coef_loud * (out_filtr(i_sample, :) .^ 2)) + ((1 - coef_loud) * mem_L_out_filtr);
    L_out_filtr(i_sample, :) = mem_L_out_filtr;
    
    % Get loudness of gain-applied outputs
    mem_L_out_L_sma = (coef_loud * (out_L_sma(i_sample, :) .^ 2)) + ((1 - coef_loud) * mem_L_out_L_sma);
    L_out_L_sma(i_sample, :) = mem_L_out_L_sma;
    mem_L_out_L_ema = (coef_loud * (out_L_ema(i_sample, :) .^ 2)) + ((1 - coef_loud) * mem_L_out_L_ema);
    L_out_L_ema(i_sample, :) = mem_L_out_L_ema;
end

%% Converts loudness obtained in log form
L_in_origin = 10 * log10(L_in_origin(:, 1));
L_out_filtr = 10 * log10(L_out_filtr(:, 1));
L_out_L_sma = 10 * log10(L_out_L_sma(:, 1));
L_out_L_ema = 10 * log10(L_out_L_ema(:, 1));

%% Plot the logarithmic (LKFS) loudness obtained

figure;
hold on;
    plot(t(1:length(t)), L_in_origin(:, 1), "LineWidth", 0.5, "Color", dict_color("yellow"))
    plot(t(1:length(t)), L_out_filtr(:, 1), "LineWidth", 0.5, "Color", dict_color("blue"))
    plot(t(1:length(t)), L_out_L_sma(:, 1), "LineWidth", 1, "Color", dict_color("green"))
    plot(t(1:length(t)), L_out_L_ema(:, 1), "LineWidth", 1, "Color", dict_color("orange"))
hold off;
ylabel("Equivalent continuous sound level (LKFS)");
xlabel("Duration /seconds");
grid on;
title("Time-Weighted L_{\tau,K} Comparisons (f_c = 440 Hz)")
ylim([-80, 00]); xlim([0, inf]);

legend(["Input (impulse)","VCF output", "SMA make-up", "EMA make-up", ], 'Location', 'northeast')
