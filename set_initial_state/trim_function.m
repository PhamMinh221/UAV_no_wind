function [state_trim, input_trim] = trim_function(UAV, Va, gamma, Radius, number_of_iteration, time_sampling)
    [alpha_trim, beta_trim, phi_trim] = minimize_J(UAV , UAV.alpha, UAV.beta, UAV.phi, Va, gamma, Radius, number_of_iteration, time_sampling),

      u_trim          = Va * cos(alpha_trim) * cos(beta_trim) ;
      v_trim          = Va * sin(beta_trim);
      w_trim         = Va * sin(alpha_trim) * cos(beta_trim);
      theta_trim    = alpha_trim + gamma ;
      psi_trim      = 0 + Va/Radius * number_of_iteration * time_sampling;
      p_trim          = - Va * sin(theta_trim) / Radius;
      q_trim          = Va/Radius * sin(phi_trim) * cos(theta_trim);
      r_trim          = Va/Radius * cos(phi_trim) * cos(theta_trim);

      state_trim = [u_trim; v_trim; w_trim; theta_trim; p_trim; q_trim; r_trim];

            %Gamma Coef.
                    GAMMA       = UAV.J_x * UAV.J_z - UAV.J_xz^2;
                    GAMMA1     = UAV.J_xz * (UAV.J_x - UAV.J_y + UAV.J_z) / GAMMA;
                    GAMMA2     = (UAV.J_z * (UAV.J_z - UAV.J_y) + UAV.J_xz ^2) / GAMMA;
                    GAMMA3     = UAV.J_z/ GAMMA;
                    GAMMA4     = UAV.J_xz / GAMMA;
                    GAMMA5     = (UAV.J_z - UAV.J_x) / UAV.J_y;
                    GAMMA6     = UAV.J_xz / UAV.J_y;
                    GAMMA7     = ((UAV.J_x - UAV.J_y)*UAV.J_x + UAV.J_xz^2)/GAMMA;
                    GAMMA8     = UAV.J_x / GAMMA;
        
                 %Dynamic Coef.
                    Cp0            = GAMMA3 * UAV.C_l0 + GAMMA4 * UAV.C_n0;
                    Cpbeta        = GAMMA3 * UAV.C_lbeta + GAMMA4 * UAV.C_nbeta;
                    Cpp            = GAMMA3 * UAV.C_lp + GAMMA4 * UAV.C_np;
                    Cpr            = GAMMA3 * UAV.C_lr + GAMMA4 * UAV.C_nr;
                    Cpdeltaa     = GAMMA3 * UAV.C_ldeltaa + GAMMA4 * UAV.C_ndeltaa;
                    Cpdeltar     = GAMMA3 * UAV.C_ldeltar + GAMMA4 * UAV.C_ndeltar;
                    Cr0            = GAMMA4 * UAV.C_l0 + GAMMA8 * UAV.C_n0;
                    Crbeta        = GAMMA4 * UAV.C_lbeta + GAMMA8 * UAV.C_nbeta;
                    Crp            = GAMMA4 * UAV.C_lp + GAMMA8 * UAV.C_np;
                    Crr            = GAMMA4 * UAV.C_lr + GAMMA8 * UAV.C_nr;
                    Crdeltaa     = GAMMA4 * UAV.C_ldeltaa + GAMMA8 * UAV.C_ndeltaa;
                    Crdeltar     = GAMMA4 * UAV.C_ldeltar + GAMMA8 * UAV.C_ndeltar;

     % trimmed input
        C_Lalpha    = UAV.C_L0 + UAV.C_Lalpha * alpha_trim ;
        C_Dalpha    = UAV.C_D0 + UAV.C_Dalpha * alpha_trim ;
        C_Xalpha    = - C_Dalpha * cos(alpha_trim) + C_Lalpha * sin(alpha_trim);
        C_Xqalpha  = -UAV.C_Dq * cos(alpha_trim) + UAV.C_Lq * sin(alpha_trim);
        C_Xdeltae   = - UAV.C_Ddeltae * cos(alpha_trim) + UAV.C_Ldeltae * sin(alpha_trim);




        elevator_trim   = ((UAV.J_xz * (p_trim^2 - r_trim)^2 + (UAV.J_x - UAV.J_z) * p_trim * r_trim)/(1/2 * UAV.rho * Va ^2 * UAV.c * UAV.S) - UAV.C_m0 - UAV.C_malpha * alpha_trim - UAV.C_mq * UAV.c * q_trim / (2*Va))/(UAV.C_mdeltae);
        thrutter_trim    = sqrt((2*UAV.m *(- r_trim * v_trim + q_trim * w_trim + UAV.gravity * sin(theta_trim)) - UAV.rho * Va^2 * UAV.S * (C_Xalpha + C_Xqalpha * UAV.c *q_trim/(2*Va) + C_Xdeltae * elevator_trim)) / (UAV.rho * UAV.S_prop * UAV.C_prop * UAV.k_motor ^2) + Va^2/(UAV.k_motor^2));
        
        sol = [Cpdeltaa, Cpdeltar; Crdeltaa, Crdeltar]^(-1) * [(-GAMMA1 * p_trim * q_trim + GAMMA2 * q_trim * r_trim)/(1/2 * UAV.rho * Va^2 * UAV.S * UAV.b) - Cp0 - Cpbeta * beta_trim - Cpp * UAV.b * p_trim /(2* Va) - Cpr * UAV.b * r_trim /(2* Va); (-GAMMA7 * p_trim * q_trim + GAMMA1 * q_trim * r_trim)/(1/2 * UAV.rho * Va^2 * UAV.S * UAV.b) - Cr0 - Crbeta * beta_trim - Crp * UAV.b * p_trim /(2* Va) - Crr * UAV.b * r_trim /(2* Va)  ]; 
        aileron_trim = sol(1);
        rudder_trim = sol(2);
        input_trim = [elevator_trim; thrutter_trim; aileron_trim; rudder_trim];
end

function [alpha, beta, phi] = minimize_J(UAV , alpha, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling)
    eps = 0.01;
    for k = 1 : 100
        alpha_plus = alpha + eps;
        beta_plus  = beta + eps ;
        phi_plus = phi +eps;
        dot_J_alpha = (compute(UAV , alpha_plus, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling) - compute(UAV , alpha, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling))/eps ;
        dot_J_beta = (compute(UAV , alpha, beta_plus, phi, Va, gamma, Radius, number_of_iteration, time_sampling) - compute(UAV , alpha, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling))/eps ;
        dot_J_phi    = (compute(UAV , alpha, beta, phi_plus, Va, gamma, Radius, number_of_iteration, time_sampling) - compute(UAV , alpha, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling))/eps ;
        alpha = alpha - 0.1 * dot_J_alpha;
        beta   = beta - 0.1 * dot_J_beta;
        phi     = phi - 0.1 * dot_J_phi;

    end


end
function J = compute(UAV , alpha, beta, phi, Va, gamma, Radius, number_of_iteration, time_sampling)
    

    dot_h_star          = Va * sin(gamma);
    dot_u_star          = 0;
    dot_v_star          = 0;
    dot_w_star          = 0;
    dot_phi_star        = 0;
    dot_theta_star     = 0;
    dot_psi_star        = Va / Radius ;
    dot_p_star          = 0;
    dot_r_star          = 0; 
    dot_q_star          = 0;
    
    dot_x_star =[dot_h_star; dot_u_star; dot_v_star; dot_w_star; dot_phi_star; dot_theta_star; dot_psi_star; dot_p_star; dot_q_star; dot_r_star];

         %Gamma Coef.
                    GAMMA       = UAV.J_x * UAV.J_z - UAV.J_xz^2;
                    GAMMA1     = UAV.J_xz * (UAV.J_x - UAV.J_y + UAV.J_z) / GAMMA;
                    GAMMA2     = (UAV.J_z * (UAV.J_z - UAV.J_y) + UAV.J_xz ^2) / GAMMA;
                    GAMMA3     = UAV.J_z/ GAMMA;
                    GAMMA4     = UAV.J_xz / GAMMA;
                    GAMMA5     = (UAV.J_z - UAV.J_x) / UAV.J_y;
                    GAMMA6     = UAV.J_xz / UAV.J_y;
                    GAMMA7     = ((UAV.J_x - UAV.J_y)*UAV.J_x + UAV.J_xz^2)/GAMMA;
                    GAMMA8     = UAV.J_x / GAMMA;
        
                 %Dynamic Coef.
                    Cp0            = GAMMA3 * UAV.C_l0 + GAMMA4 * UAV.C_n0;
                    Cpbeta        = GAMMA3 * UAV.C_lbeta + GAMMA4 * UAV.C_nbeta;
                    Cpp            = GAMMA3 * UAV.C_lp + GAMMA4 * UAV.C_np;
                    Cpr            = GAMMA3 * UAV.C_lr + GAMMA4 * UAV.C_nr;
                    Cpdeltaa     = GAMMA3 * UAV.C_ldeltaa + GAMMA4 * UAV.C_ndeltaa;
                    Cpdeltar     = GAMMA3 * UAV.C_ldeltar + GAMMA4 * UAV.C_ndeltar;
                    Cr0            = GAMMA4 * UAV.C_l0 + GAMMA8 * UAV.C_n0;
                    Crbeta        = GAMMA4 * UAV.C_lbeta + GAMMA8 * UAV.C_nbeta;
                    Crp            = GAMMA4 * UAV.C_lp + GAMMA8 * UAV.C_np;
                    Crr            = GAMMA4 * UAV.C_lr + GAMMA8 * UAV.C_nr;
                    Crdeltaa     = GAMMA4 * UAV.C_ldeltaa + GAMMA8 * UAV.C_ndeltaa;
                    Crdeltar     = GAMMA4 * UAV.C_ldeltar + GAMMA8 * UAV.C_ndeltar;
    % trimmed state
     
        UAV.u          = Va * cos(alpha) * cos(beta) ;
        UAV.v          = Va * sin(beta);
        UAV.w         = Va * sin(alpha) * cos(beta);
        UAV.theta    = alpha + gamma ;
        UAV.psi        = 0 + Va/Radius * number_of_iteration * time_sampling;
        UAV.p          = - Va * sin(UAV.theta) / Radius;
        UAV.q          = Va/Radius * sin(phi) * cos(UAV.theta);
        UAV.r          = Va/Radius * cos(phi) * cos(UAV.theta);
        
     % trimmed input
        C_Lalpha    = UAV.C_L0 + UAV.C_Lalpha * alpha ;
        C_Dalpha    = UAV.C_D0 + UAV.C_Dalpha * alpha ;
        C_Xalpha    = - C_Dalpha * cos(alpha) + C_Lalpha * sin(alpha);
        C_Xqalpha  = -UAV.C_Dq * cos(alpha) + UAV.C_Lq * sin(alpha);
        C_Xdeltae   = - UAV.C_Ddeltae * cos(alpha) + UAV.C_Ldeltae * sin(alpha);




        elevator_star   = ((UAV.J_xz * (UAV.p^2 - UAV.r)^2 + (UAV.J_x - UAV.J_z) * UAV.p * UAV.r)/(1/2 * UAV.rho * Va ^2 * UAV.c * UAV.S) - UAV.C_m0 - UAV.C_malpha * alpha - UAV.C_mq * UAV.c * UAV.q / (2*Va))/(UAV.C_mdeltae);
        thrutter_star    = sqrt((2*UAV.m *(- UAV.r * UAV.v + UAV.q * UAV.w + UAV.gravity * sin(UAV.theta)) - UAV.rho * Va^2 * UAV.S * (C_Xalpha + C_Xqalpha * UAV.c *UAV.q/(2*Va) + C_Xdeltae * elevator_star)) / (UAV.rho * UAV.S_prop * UAV.C_prop * UAV.k_motor ^2) + Va^2/(UAV.k_motor^2));
        
        sol = [Cpdeltaa, Cpdeltar; Crdeltaa, Crdeltar]^(-1) * [(-GAMMA1 * UAV.p * UAV.q + GAMMA2 * UAV.q * UAV.r)/(1/2 * UAV.rho * Va^2 * UAV.S * UAV.b) - Cp0 - Cpbeta * beta - Cpp * UAV.b * UAV.p /(2* Va) - Cpr * UAV.b * UAV.r /(2* Va); (-GAMMA7 * UAV.p * UAV.q + GAMMA1 * UAV.q * UAV.r)/(1/2 * UAV.rho * Va^2 * UAV.S * UAV.b) - Cr0 - Crbeta * beta - Crp * UAV.b * UAV.p /(2* Va) - Crr * UAV.b * UAV.r /(2* Va)  ]; 
        aileron_star = sol(1);
        rudder_star = sol(2);
        
        f = UAV.find_f(elevator_star, thrutter_star, aileron_star, rudder_star);

        J = norm(dot_x_star - f)^2;

end 