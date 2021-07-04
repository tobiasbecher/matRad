function [DoseRate] = getDoseRate2D(obj,r,theta)
% This function returns the dose rate according to the 2D formalism 
% in Rivard et al. (2004): AAPM TG-43 update, page 637, Eq. 1.
%
% Output parameters:
% DoseRate: Dose rate in water at P(r, tetha). dose rate is returned in
%           Gy/s. 
%           When a time is specified, Dose rate is returned in Gy i.e. the
%           dose is returned !!!
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to 
%           P(r,theta), with units of cm.
% theta:    column vector of polar angles between the longitudinal axis of the source 
%           and the ray from the active source center to the calculation
%           point, units of 
% time:     treatment time in s

global cm
cm = 1;
% TODO
% r0: The reference distance, which is 1 cm for this protocol
r0       = 1*cm;

% Sk: Air-kerma strength
Sk       = obj.SourceStrengthImplanted;
% lambda: Dose-rate constant in water (Lambda)
lambda   = obj.lambda;
% GL0: Geometry function G_L(r_0, theta_0)
GL0      = obj.getTransverseGeometryFunction(r0);
% GL: Geometry function G_L(r, theta) 
GL       = obj.getGeometryFunction(r,theta);
% F: 2D anisotropy function F(r,theta)
F        = obj.get2DAnisotropyFunction(r,theta);
% gl: radial dose function g_L(r)
gl       = obj.getRadialDoseFunction(r);


% 2D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
% page 637, Eq. 1:
DoseRate = Sk*lambda*((GL./GL0).*gl.*F);


end

