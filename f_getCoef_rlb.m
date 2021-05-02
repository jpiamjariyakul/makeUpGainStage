%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revised low-frequency B-curve filter coefficients
%
% Author: Jay Piamjariyakul
%
% Sources
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