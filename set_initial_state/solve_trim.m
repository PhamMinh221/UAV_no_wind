function [state_trim , input_trim] = solve_trim(self, Va_command, R_command)
    
        alpha = 0; beta = 0; phi = 0; eps = 0.0001;
        for k = 1 : 2500
            alpha_plus = alpha + eps;
            beta_plus  = beta + eps;
            phi_plus    = phi + eps;

            dot_J_alpha = (compute_J(self, alpha_plus, beta, phi, Va_command, R_command) - compute_J(self, alpha, beta, phi, Va_command, R_command)) ;
            dot_J_beta   = (compute_J(self, alpha, beta_plus, phi, Va_command, R_command) - compute_J(self, alpha, beta, phi, Va_command, R_command)) ;
            dot_J_phi     = (compute_J(self, alpha, beta, phi_plus, Va_command, R_command) - compute_J(self, alpha, beta, phi, Va_command, R_command));

            Gd = 0.004;
            alpha = alpha - Gd * dot_J_alpha;
            beta  = beta - Gd* dot_J_beta;
            phi    = phi - Gd * dot_J_phi   ;
        end

        alpha_trim          = alpha;
        beta_trim           = beta;
        phi_trim             = phi;
        u_trim                =  Va_command * cos(alpha) * cos (beta);
        v_trim                = Va_command * sin(beta);
        w_trim               = Va_command * sin(alpha) * cos(beta);
        theta_trim          = alpha;
        p_trim                = -Va_command/R_command * sin(theta_trim);
        q_trim                = Va_command/R_command * sin(phi) * cos(theta_trim);
        r_trim                 = Va_command/R_command * cos(phi) * cos(theta_trim);

       state_trim  = [u_trim, v_trim, w_trim, theta_trim, p_trim, q_trim, r_trim, alpha_trim, beta_trim, phi_trim];


        delta_e        = ((self.J_xz*(p_trim^2 - r_trim^2)+ (self.J_x -self.J_z) * p_trim * r_trim) / (0.5 * self.rho * Va_command^2 * self.c * self.S) - self.C_m0 - self.C_malpha * alpha_trim - self.C_mq * self.c * q_trim / (2 * Va_command))/ self.C_mdeltae;
        
        CLalpha     = self.C_L0 + self.C_Lalpha * alpha_trim;
        CDalpha     = self.C_D0 + self.C_Dalpha *alpha_trim;
        CXalpha     = - CDalpha * cos(alpha_trim) + CLalpha * sin(alpha_trim);
        CXqalpha   = - self.C_Dq * cos(alpha_trim) + self.C_Lq * sin(alpha_trim);
        CXdeltae    = - self.C_Ddeltae * cos(alpha_trim) + self.C_Ldeltae * sin(alpha_trim);
  
        delta_t      = sqrt((2 * self.m * (-r_trim * v_trim + q_trim * w_trim + self.gravity * sin(theta_trim)) - self.rho * Va_command^2 * self.S * (CXalpha + CXqalpha * self.c * q_trim /(2* Va_command) + CXdeltae* delta_e))/(self.rho* self.S_prop * self.C_prop * self.k_motor^2) + Va_command^2/self.k_motor^2);

        sol = [self.C_pdeltaa, self.C_pdeltar; self.C_rdeltaa, self.C_rdeltar]^(-1) * [( - self.GAMMA1 * p_trim * q_trim + self.GAMMA2 * q_trim * r_trim )/(0.5 *self.rho * Va_command * Va_command * self.S * self.b ) - self.C_p0 - self.C_pbeta * beta_trim - self.C_pp * self.b * self.p/(2 * Va_command) - self.C_pr * self.b * r_trim /(2* Va_command); (-self.GAMMA7 * p_trim * q_trim + self.GAMMA1 * q_trim * r_trim)/(0.5 *self.rho * Va_command^2 * self.S * self.b) - self.C_r0 - self.C_rbeta * beta_trim - self.C_rp * self.b * p_trim /(2 * Va_command) - self.C_rr * self.b * r_trim /(2 * Va_command) ];

        delta_a = sol(1);
        delta_r = sol(2);

        input_trim = [delta_e, delta_t, delta_a, delta_r];

end

function J = compute_J(self, alpha, beta, phi, Va_command, R_command)
    
        dot_x_trim =[0; 0; 0; 0; 0; 0; Va_command/R_command; 0; 0; 0];
        
        self.alpha     = alpha;
        self.phi       = phi;
        self.beta      = beta;
        self.u          = Va_command * cos(alpha) * cos (beta);
        self.v          = Va_command * sin(beta);
        self.w         = Va_command * sin(alpha) * cos(beta);
        self.theta    = alpha;
        self.p          = -Va_command/R_command * sin(self.theta);
        self.q          = Va_command/R_command * sin(self.phi) * cos(self.theta);
        self.r           = Va_command/R_command * cos(self.phi) * cos(self.theta);

        delta_e        = ((self.J_xz*(self.p^2 - self.r^2)+ (self.J_x -self.J_z) * self.p * self.r) / (0.5 * self.rho * Va_command^2 * self.c * self.S) - self.C_m0 - self.C_malpha * self.alpha - self.C_mq * self.c *self.q / (2 * Va_command))/ self.C_mdeltae;
        
        CLalpha     = self.C_L0 + self.C_Lalpha * self.alpha;
        CDalpha     = self.C_D0 + self.C_Dalpha *self.alpha;
        CXalpha     = - CDalpha * cos(self.alpha) + CLalpha * sin(self.alpha);
        CXqalpha   = - self.C_Dq * cos(self.alpha) + self.C_Lq * sin(self.alpha);
        CXdeltae    = - self.C_Ddeltae * cos(self.alpha) + self.C_Ldeltae * sin(self.alpha);
  
        delta_t      = sqrt((2 * self.m * (-self.r * self.v + self.q * self.w + self.gravity * sin(self.theta)) - self.rho * Va_command^2 * self.S * (CXalpha + CXqalpha * self.c * self.q /(2* Va_command) + CXdeltae* delta_e))/(self.rho* self.S_prop * self.C_prop * self.k_motor^2) + Va_command^2/self.k_motor^2);

        sol = [self.C_pdeltaa, self.C_pdeltar; self.C_rdeltaa, self.C_rdeltar]^(-1) * [( - self.GAMMA1 * self.p * self.q + self.GAMMA2 * self.q * self.r )/(0.5 *self.rho * Va_command* self.S * self.b ) - self.C_p0 - self.C_pbeta * self.beta - self.C_pp * self.b * self.p/(2 * Va_command) - self.C_pr * self.b * self.r /(2* Va_command); (-self.GAMMA7 * self.p * self.q + self.GAMMA1 * self.q * self.r)/(0.5 *self.rho * Va_command^2 * self.S * self.b) - self.C_r0 - self.C_rbeta * self.beta - self.C_rp * self.b * self.p /(2 * Va_command) - self.C_rr * self.b * self.r /(2 * Va_command) ];

        delta_a = sol(1);
        delta_r = sol(2);

        f = self.find_f(delta_e, delta_t, delta_a, delta_r);

        J = norm(dot_x_trim - f)^2;

end