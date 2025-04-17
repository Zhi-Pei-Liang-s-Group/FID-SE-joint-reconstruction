function sxt_recon_final = Jresolved_MRSI_FIDSE_jointRecon(sxtc_meta_multiTE,Vt_meta,mu_coef,var_coef,fm_hz_wat,sense_map_odd,tvec_recon_cell,brainMask,anatRef)
% This is the packaged code for joint FID/SE-reconstruction for J-resolved MRSI
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
%   sxt_recon_final:    final reconstruction
%
% 
%   Liang's Lab @ UIUC
%   Created: 2022-12-20
%

    %% ------------ step #01: parse the input --------------- %%

    %% image size
    N_TE              = length(sxtc_meta_multiTE);
    [L1,L2,L3]        = size_dims(sxtc_meta_multiTE{1},1:3);
    M_cell            = cell(N_TE,1);
    for index_te = 1:N_TE
        M_cell{index_te} = size(sxtc_meta_multiTE{index_te},4);
    end

    %% ------------ step #02: derive lambda based on distribution --------------- %%

    %% coil combination
    sxt_meta_multiTE  = cell(N_TE,1);
    sxt_meta_multiTE  = sxtc_meta_multiTE;
    for index_te = 2:N_TE
        sxt_tmp = sxt_meta_multiTE{index_te};
        sxt_tmp = combineData_SenseMap(imresize_ktrunc(sxt_tmp,size_dims(sense_map_odd,1:3),0,0),sense_map_odd);
        sxt_meta_multiTE{index_te} = sxt_tmp;
    end

    %% noise estimation
    sxt_spice_low   = imresize_ktrunc(sxt_meta_multiTE{1},[36,36,12],false);
    noise_samp      = vec(sxt_spice_low(:,:,:,end-10:end));
    noiseVar        = var(noise_samp);

    %% MAP
    % noise level
    noise_var       = noiseVar;

    % distribution parameter
    mu_meta         = vec(mu_coef);
    lambda_meta     = vec(noise_var./var_coef);

    %% ---------- Step #3: MAP-based recon 1st iteration ----------%%

    %% prepare the data and sampling masks
    dkt_cell           = cell(N_TE,1);
    ktMask_cell        = cell(N_TE,1);
    for index_te = 1:N_TE
        dkt_cell{index_te}    = F3_x2k(sxtc_meta_multiTE{index_te});
        ktMask_cell{index_te} = ktrunc(true(size_dims(sxtc_meta_multiTE{index_te},1:4)),[L1,L2,L3]);
    end

    %% anatomical weights
    Iref                = normRange(imresize3d(anatRef,[L1,L2,L3]));
    opt                 = struct('mask',brainMask,'verbose',0);
    Dg                  = genWL2mtx(Iref,opt);

    %% recon parameters
    lambda_wl2          = 0.5; % selected by discrepancy principle
    ind_fit             = 1:sum(cell2mat(M_cell));

    %% reconstruction reconstruction
    opt = struct('lambda_wl2',lambda_wl2,'ind_fit',ind_fit,'MaxIter',20,'disp_flag',false,'sxt_env',[]);
    sxt_recon_meta = SPICErecon_multiCoil_multiTE_wGauss_wEnv_wl2(dkt_cell,ktMask_cell,fm_hz_wat,sense_map_odd,brainMask,tvec_recon_cell,Vt_meta,mu_meta,lambda_meta*0.1,Dg,opt);
    sxt_recon = sxt_recon_meta;

    %% -----------  Step #4: recon the residual --------- %%

    %% residual data
    dkt_res_cell              = cell(N_TE,1);
    Ind                       = 1;
    for index_te = 1:N_TE
        dkt_tmp               = dkt_cell{index_te};
        tmp_recon             = sxt_recon_meta(:,:,:,Ind:Ind+M_cell{index_te}-1);
        tmp_recon             = ConjugatePhase3D(tmp_recon,-fm_hz_wat,tvec_recon_cell{index_te});
        if index_te == 1
            tmp_recon         = F3_x2k(tmp_recon);
            tmp_recon         = oper4D_mask_4D_data(tmp_recon,ktMask_cell{index_te}); % sampling mask
        else
            tmp_recon         = bsxfun(@times,tmp_recon,sense_map_odd); % sensitivity map
            tmp_recon         = F3_x2k(tmp_recon);
            tmp_recon         = oper4D_mask_4D_data(tmp_recon,ktMask_cell{index_te}); % sampling mask
        end
        dkt_res_cell{index_te} = dkt_tmp - reshape(tmp_recon,size(dkt_tmp));
        Ind                   = Ind+M_cell{index_te};
    end

    %% derive the residual subspace
    low_res          = [10,10,6];
    sxt_res          = [];
    for index_te = 1:N_TE
        if index_te == 1
            tmp = F3_k2x(dkt_res_cell{index_te});
            tmp = ConjugatePhase3D(tmp,fm_hz_wat,tvec_recon_cell{index_te});
        else
            tmp = F3_k2x(ktrunc(dkt_res_cell{index_te},[L1,L2,L3]));
            tmp = combineData_SenseMap(tmp,sense_map_odd);
            tmp = ConjugatePhase3D(tmp,fm_hz_wat,tvec_recon_cell{index_te});
        end
        sxt_res = cat(4,sxt_res,imresize_ktrunc(tmp,low_res,0,0));
    end

    %% residual subspace
    [U_res,S_res,V_res] = estMathBases(sxt_res,[]);

    %% add weighting
    Rank_res          = 8; % selected based on the residual singular values
    bases_weight      = [((diag(S_res(1:Rank_res,1:Rank_res))/S_res(1,1))).'];
    V_res_weight      = bsxfun(@times,[V_res(:,1:Rank_res)],bases_weight);

    %% residual reconstruction
    lambda_wl2        = 3;
    opt = struct('lambda_wl2',lambda_wl2,'ind_fit',ind_fit,'MaxIter',20,'disp_flag',false,'sxt_env',[]);
    [sxt_recon_res]   = SPICErecon_multiCoil_multiTE_wGauss_wEnv_wl2(dkt_res_cell,ktMask_cell,fm_hz_wat,sense_map_odd,brainMask,tvec_recon_cell,V_res_weight,zeros(Rank_res,1,'single'),zeros(Rank_res,1,'single'),Dg,opt);


    %% -----------  Step #5: organize the outputs --------- %%

    %% combine recon results
    sxt_recon    = sxt_recon + sxt_recon_res;
    sxt_recon_final  = cell(N_TE,1);
    Ind = 1;
    for index_te = 1:N_TE
        sxt_recon_final{index_te}  = sxt_recon(:,:,:,Ind:Ind+M_cell{index_te}-1);
        Ind = Ind+M_cell{index_te};
    end

end
