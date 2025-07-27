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
  
        delta_t      = sqrt(2 * self.m * (-self.r * self.v + self.q * self.w + self.gravity * sin(self.theta)) - self.rho * Va_command^2 * self.S * (CXalpha + CXqalpha * self.c * self.q /(2* Va_command) + CXdeltae* delta_e)/(self.rho* self.S_prop * self.C_prop * self.k_motor^2) + Va_command^2/self.k_motor^2);

        sol = [self.C_pdeltaa, self.C_pdeltar; self.C_rdeltaa, self.C_rdeltar]^(-1) * [( - self.GAMMA1 * self.p * self.q + self.GAMMA2 * self.q * self.r )/(0.5 *self.rho * Va_command* self.S * self.b ) - self.C_p0 - self.C_pbeta * self.beta - self.C_pp * self.b * self.p/(2 * Va_command) - self.C_pr * self.b * self.r /(2* Va_command); (-self.GAMMA7 * self.p * self.q + self.GAMMA1 * self.q * self.r)/(0.5 *self.rho * Va_command^2 * self.S * self.b) - self.C_r0 - self.C_rbeta * self.beta - self.C_rp * self.b * self.p /(2 * Va_command) - self.C_rr * self.b * self.r /(2 * Va_command) ];

        delta_a = sol(1);
        delta_r = sol(2);

        f = self.find_f(delta_e, delta_t, delta_a, delta_r);

        J = norm(dot_x_trim - f)^2;

end