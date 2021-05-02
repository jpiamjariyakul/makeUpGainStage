%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-Dimensional Biquad filter
%
% Author: Jay Piamjariyakul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_filt, seq_delay] = f_1dFilter(b, a, in_audio, seq_delay)

b0 = b(1);
b1 = b(2);
b2 = b(3);
a1 = a(2);
a2 = a(3);

sz_channel = size(in_audio, 2);
out_filt = zeros(size(in_audio, 1), sz_channel);

for i_sample = 1:size(in_audio, 1)
delayNow = in_audio(i_sample, :) - (a1 * seq_delay(1, :)) - (a2 * seq_delay(2, :));

out_filt(i_sample, :) = (b0 * delayNow) + (b1 * seq_delay(1, :)) + (b2 * seq_delay(2, :));

seq_delay(2, :) = seq_delay(1, :);
seq_delay(1, :) = delayNow;

end

end