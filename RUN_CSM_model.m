%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% RUN script for Comprehensive SNc model


% SNc with ATP dynamics (Francis et.al., 2013)
% Dopamine synthesis, storage, release, metabolism and terminal autoreceptors (Bravo, 2012)
% Ca2+ induced apoptosis (Hong et.al., 2012)
% Calcium-induced calcium release (Marhl et.al., 2000)
% Energy Metabolism (Cloutier & Wellstead, 2010)
% PD pathology pathways (Cloutier & Wellstead, 2012)

%% CODE

gl = 1; % Glucose concentration in mM
mt = 1; % Extend of oxygen available (0-no oxygen; 1-adequate oxygen)
filename = 'test';

dur=1000; % Duration of simulation in milliseconds

VTA_ATPapopNM(dur,gl,mt,filename)