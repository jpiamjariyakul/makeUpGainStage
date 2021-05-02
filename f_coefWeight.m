%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential moving average filter - weighting coefficient
%
% Author: Jay Piamjariyakul
%
% Sources
% EMA weighting coefficient:
% - D. Ward, "Applications of loudness models in audio engineering," Ph.D. 
%   dissertation, Birmingham City University, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coefWeight = f_coefWeight(Fc, Fs)
coefWeight = 1 - exp( (-2 * pi) * ( Fc / Fs ) );
end