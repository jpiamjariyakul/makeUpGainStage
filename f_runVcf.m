%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input buffer generator for f_vcf_nonlinear.m
%
% Author: Jay Piamjariyakul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = f_runVcf(in_audio, g, k)
y = zeros(4, size(in_audio, 2));
y_d2 = 0;

out = zeros(size(in_audio, 1), 2);

for i_sample = 1:size(in_audio, 1)
    [y, y_d2] = f_vcf_nonlinear(y, in_audio(i_sample, :), g, k, y_d2);
    out(i_sample, :) = y(4, :);
end

end