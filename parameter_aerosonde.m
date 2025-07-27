
%% Parameters for the Aerosonde UAV
    MAV.m            = 1.56;                 %[kg]
    MAV.J_x          = 0.1147;              %[kg-m^2]
    MAV.J_y          = 0.0576;               %[kg-m^2]
    MAV.J_z          = 0.1712;               %[kg-m^2]
    MAV.J_xz        = 0.0015;                %[kg-m^2]
    MAV.S            = 0.2589;                %[m^2]
    MAV.b            = 1.4224;                  %[m]
    MAV.c             = 0.3302;                 %[m]
    MAV.S_prop    = 0.0314;                 %[m^2]
    MAV.rho          = 1.2682;                 %[kg/m^3]
    MAV.k_motor   = 80;                        
    MAV.k_Tp        = 0;
    MAV.k_w         = 0;
    MAV.gravity     =9.81;

 %% Aerodynamic coefficients for the Aerosonde UAV
    % Longitudinal Coef.
    MAV.C_L0                = 0.28;
    MAV.C_D0                = 0.03;
    MAV.C_m0                = -0.02338;
    MAV.C_Lalpha            = 3.45;
    MAV.C_Dalpha            = 0.30;
    MAV.C_malpha            = -0.38;
    MAV.C_Lq                   = 0;
    MAV.C_Dq                   = 0;
    MAV.C_mq                   = -3.6;
    MAV.C_Ldeltae             = -0.36;
    MAV.C_Ddeltae             = 0;
    MAV.C_mdeltae            = -0.5;
    MAV.C_prop                 = 1.0 ;
    MAV.M                         = 50 ;
    MAV.alpha0                  = 0.4712;
    MAV.epsilon                 = 0.1592;
    % Lateral Coef.
    MAV.C_Y0                    = 0 ;
    MAV.C_l0                     = 0 ;
    MAV.C_n0                    = 0 ;
    MAV.C_Ybeta                = -0.98 ;
    MAV.C_lbeta                 = -0.12 ;
    MAV.C_nbeta                = 0.25 ;
    MAV.C_Yp                     = 0 ;
    MAV.C_lp                      = -0.26 ;
    MAV.C_np                     = 0.022 ;
    MAV.C_Yr                      = 0;
    MAV.C_lr                       = 0.14;
    MAV.C_nr                      = -0.35;
    MAV.C_Ydeltaa               = 0;
    MAV.C_ldeltaa               = 0.08;
    MAV.C_ndeltaa              = 0.06;
    MAV.C_Ydeltar             = -0.17;
    MAV.C_ldeltar               = 0.105;
    MAV.C_ndeltar              = -0.032;

    %  GAMMA Coef.
    MAV.GAMMA       = MAV.J_x * MAV.J_z - MAV.J_xz^2;
    MAV.GAMMA1     = MAV.J_xz * (MAV.J_x - MAV.J_y + MAV.J_z) / MAV.GAMMA;
    MAV.GAMMA2     = (MAV.J_z * (MAV.J_z - MAV.J_y) + MAV.J_xz ^2) / MAV.GAMMA;
    MAV.GAMMA3     = MAV.J_z/ MAV.GAMMA;
    MAV.GAMMA4     = MAV.J_xz / MAV.GAMMA;
    MAV.GAMMA5     = (MAV.J_z - MAV.J_x) / MAV.J_y;
    MAV.GAMMA6     = MAV.J_xz / MAV.J_y;
    MAV.GAMMA7     = ((MAV.J_x - MAV.J_y)*MAV.J_x + MAV.J_xz^2)/MAV.GAMMA;
    MAV.GAMMA8     = MAV.J_x / MAV.GAMMA;

    %
    MAV.C_p0                     = MAV.GAMMA3 * MAV.C_l0 + MAV.GAMMA4 * MAV.C_n0;
    MAV.C_pbeta                 = MAV.GAMMA3 * MAV.C_lbeta + MAV.GAMMA4 * MAV.C_nbeta;
    MAV.C_pp                      = MAV.GAMMA3 * MAV.C_lp + MAV.GAMMA4 * MAV.C_np;
    MAV.C_pr                      = MAV.GAMMA3 * MAV.C_lr + MAV.GAMMA4 * MAV.C_nr;
    MAV.C_pdeltaa               = MAV.GAMMA3 * MAV.C_ldeltaa + MAV.GAMMA4 * MAV.C_ndeltaa;
    MAV.C_pdeltar                = MAV.GAMMA3 * MAV.C_ldeltar + MAV.GAMMA4 * MAV.C_ndeltar;
    MAV.C_r0                       = MAV.GAMMA4 * MAV.C_l0 + MAV.GAMMA8 * MAV.C_n0;
    MAV.C_rbeta                   = MAV.GAMMA4 * MAV.C_lbeta + MAV.GAMMA8 * MAV.C_nbeta;
    MAV.C_rp                       = MAV.GAMMA4 * MAV.C_lp + MAV.GAMMA8 * MAV.C_np;
    MAV.C_rr                        = MAV.GAMMA4 * MAV.C_lr + MAV.GAMMA8 * MAV.C_nr;
    MAV.C_rdeltaa                = MAV.GAMMA4 * MAV.C_ldeltaa + MAV.GAMMA8 * MAV.C_ndeltaa;
    MAV.C_rdeltar                 = MAV.GAMMA4 * MAV.C_ldeltar + MAV.GAMMA8 * MAV.C_ndeltar;
    
    
