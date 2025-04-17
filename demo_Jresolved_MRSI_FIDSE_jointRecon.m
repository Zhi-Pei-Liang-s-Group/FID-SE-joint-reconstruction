% This is the demonstration for joint FID/SE joint reconstruction for J-resolved MRSI
% 
% [Inputs]
%   sxtc_meta_multiTE:  cell of input data, each cell contains the nuisance removed k-space data for FID or SE
%   Vt_meta:            joint recon subspace for FID/SE
%   mu_coef:            mean of prior distribution
%   var_coef:           variance of prior distribution
%   fm_hz_wat:          B0 field map
%   sense_map_odd:      coil sensitivity map
%   tvec_recon_cell:    cell of time vectors corresponding to 'sxtc_meta_multiTE'
%   brainMask:          binary mask for brain support
%   anatRef:            anatomical reference image
%
% [Outputs]
%   sxt_recon_final:    final reconstruction results
%
% 
%   Liang's Lab @ UIUC
%   Created: 2022-12-20
%
clear;
restoredefaultpath;

addpath('./pcode');
load('./data/demo_jointRecon.mat');

sxt_recon_final = Jresolved_MRSI_FIDSE_jointRecon(sxtc_meta_multiTE,Vt_meta,mu_coef,var_coef,fm_hz_wat,sense_map_odd,tvec_recon_cell,brainMask,anatRef);
