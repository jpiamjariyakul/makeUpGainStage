%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear Moog filter - Difference equations
%
% Author: Jay Piamjariyakul
%
% Sources
% Unweighted nonlinear digital Moog VCF difference equations:
% - P. Daly, "A comparison of virtual analogue Moog VCF models," Master's
%   thesis, Univ. ofEdinburgh, Edinburgh, UK, Aug, 2012.
% Original nonlinear implementation:
% - A. Huovilainen, "Non-linear digital implementation of the Moog ladder
%   filter," in Proceed-ings of the International Conference on Digital 
%   Audio Effects (DAFx-04), 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, y_d2] = f_filter_nonlinear(y, in, g, k, y_d2)
% Nonlinear implementation at y(4)+y_d2 uses previous value stored
% Initial previous value (ie y_d2) = 0
y_0 = in - ( k * ( 0.5 * ( y(4, :) + y_d2 ) ) );
y(1, :) = y(1, :) + ( g * ( tanh( y_0  ) - tanh( y(1, :) ) ) );
y(2, :) = y(2, :) + ( g * ( tanh( y(1, :) ) - tanh( y(2, :) ) ) );
y(3, :) = y(3, :) + ( g * ( tanh( y(2, :) ) - tanh( y(3, :) ) ) );
y_d2 = y(4, :);
y(4, :) = y(4, :) + ( g * ( tanh( y(3, :) ) - tanh( y(4, :) ) ) );
end