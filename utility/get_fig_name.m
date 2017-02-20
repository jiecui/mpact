function fig_name = get_fig_name(ref_alg)
% GET_FIG_NAME (summary)
%
% Syntax:
%
% Input(s):
%
% Output(s):
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2017 Richard J. Cui. Created: Sun 02/19/2017 10:20:53.731 PM
% $Revision: 0.1 $  $Date: Sun 02/19/2017 10:20:53.739 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

switch ref_alg
    case 'ExpectMax'
        fig_name = 'Coarse: ChirpletMP + Refine: EM';
    case 'MaxLikeliEst'
        fig_name = 'Coarse: ChirpletMP + Refine: MLE';
    otherwise
        error('sim2_x_chirps:mp_act_signal', 'Unknow refinement algorithm')
    
end % switch

end % function get_fig_name

% [EOF]
