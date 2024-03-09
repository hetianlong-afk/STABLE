% Longitudinal monopole HOMs from main or harmonic cavities
% input in column
%% parameters:
% case: W/O HOMs
Q_hom_m0 = [];
R_hom_m0 = [];
fre_hom_m0 = [];

% case: One HOM
% Q_hom_m0 = [100];
% R_hom_m0 = [50]*Q_hom_m0;
% fre_hom_m0 = [299792458/480*1200+100e3];

% CESR-B HOMs
% Q_hom_m0 = [350;420;150;60;20;40;30];
% R_hom_m0 = [8.14;54.6;505;13.2;65.6;46;7.5];
% fre_hom_m0 = [950.55;976.62;1014.38;1181.5;1361;1481.5;1580]*1e6;

%% Add to HALF structure
if ~isempty(Q_hom_m0)
HALF.Q_hom_m0 = Q_hom_m0;
HALF.R_hom_m0 = R_hom_m0;
HALF.fre_hom_m0 = fre_hom_m0;

HALF.wrf_hom_m0 = 2*pi*HALF.fre_hom_m0;
HALF.angle_hom_m0 = HALF.wrf_hom_m0 * HALF.Tb;

HALF.VbImagFactor_hom_m0   = 1./sqrt(4*HALF.Q_hom_m0.^2-1);
HALF.rot_coef_hom_m0      = sqrt(1-1./(2*HALF.Q_hom_m0).^2);

HALF.Q_hom_length = length(HALF.Q_hom_m0);
end