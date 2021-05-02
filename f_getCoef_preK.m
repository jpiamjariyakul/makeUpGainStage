%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-K-weighting filter coefficients
%
% Author: Jay Piamjariyakul
%
% Sources
% Method for obtaining parameter values at any fs:
% - D. Ward, "Applications of loudness models in audio engineering," Ph.D. 
%   dissertation, Birmingham City University, 2017
% Original parameter values for fs=48000 Hz:
% - International Telecommunications Union, "ITU-R BS.1770 Algorithms to 
%   measure audioprogramme loudness and true-peak audio level," 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fil_preK_coef_b, fil_preK_coef_a] = f_getCoef_preK(Fs_filter)
fil_preK_fc = 1681.9744510;
fil_preK_q = 0.7071752;
fil_preK_VL = 1;
fil_preK_VB = 1.2587209;
fil_preK_VH = 1.5848647;
fil_preK_Omega = tan(pi * (fil_preK_fc / Fs_filter)); % The only fs-dependent param

fil_preK_coef_a = [  (fil_preK_Omega ^ 2) + (fil_preK_Omega / fil_preK_q) + 1,  ... % Coefficient a0
                    2 * ((fil_preK_Omega^2) - 1),    ...                         % Coefficient a1
                    (fil_preK_Omega^2) - (fil_preK_Omega / fil_preK_q) + 1 ...     % Coefficient a2
                    ];
fil_preK_coef_b = [  (fil_preK_VL * fil_preK_Omega^2) + (fil_preK_VB * fil_preK_Omega / fil_preK_q) + fil_preK_VH, ... % Coefficient b0
                    2 * ((fil_preK_VL * fil_preK_Omega^2) - fil_preK_VH),  ...                                     % Coefficient b1
                    (fil_preK_VL * fil_preK_Omega^2) - (fil_preK_VB * fil_preK_Omega / fil_preK_q) + fil_preK_VH  ... % Coefficient b2
                    ];
fil_preK_coef_a0 = fil_preK_coef_a(1);
fil_preK_coef_a = fil_preK_coef_a / fil_preK_coef_a0; % Normalise coefficients to a0
fil_preK_coef_b = fil_preK_coef_b / fil_preK_coef_a0; % Normalise coefficients to a0
end