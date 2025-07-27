classdef UAV 

    properties 
       m, J_x, J_y, J_z, J_xz, S, b, c, S_prop, rho, k_motor, k_Tp, k_w, gravity, Va, f, 

       % Longitudinal Coef.
       C_L0, C_D0, C_m0, C_Lalpha, C_Dalpha, C_malpha, C_Lq, C_Dq, C_mq, C_Ldeltae, C_Ddeltae, C_mdeltae, C_prop, M, alpha0, epsilon,
       % Lateral Coef.
       C_Y0, C_l0, C_n0, C_Ybeta, C_lbeta, C_nbeta, C_Yp, C_lp, C_np,
       C_Yr, C_lr, C_nr, C_Ydeltaa, C_ldeltaa, C_ndeltaa, C_Ydeltar, C_ldeltar, C_ndeltar,

       pn, pe, pd, u, v, w, phi, theta, psi, p, q, r, alpha, beta, chi, gamma, chic, gammaa

       GAMMA , GAMMA1, GAMMA2, GAMMA3, GAMMA4, GAMMA5, GAMMA6, GAMMA7, GAMMA8,
       C_p0, C_pbeta, C_pp, C_pr, C_pdeltaa, C_pdeltar, C_r0, C_rbeta, C_rp, C_rr, C_rdeltaa, C_rdeltar
    end

    methods

        %dynamic system
        function self = update_state(self, u1, u2, u3, u4, u5)

                 %Dynamic Coef.
                    CLalpha     = self.C_L0 + self.C_Lalpha * self.alpha;
                    CDalpha     = self.C_D0 + self.C_Dalpha *self.alpha;
                    CXalpha     = - CDalpha * cos(self.alpha) + CLalpha * sin(self.alpha);
                    CXqalpha   = - self.C_Dq * cos(self.alpha) + self.C_Lq * sin(self.alpha);
                    CXdeltae    = - self.C_Ddeltae * cos(self.alpha) + self.C_Ldeltae * sin(self.alpha);
                    CZalpha     = - CDalpha * sin(self.alpha) - CLalpha * cos(self.alpha);
                    CZqalpha   = - self.C_Dq * sin(self.alpha) - self.C_Lq * cos(self.alpha);
                    CZdeltae    = - self.C_Ddeltae * sin(self.alpha) - self.C_Ldeltae * cos(self.alpha);

                % deriaviate of state
                    dot_pn      = cos(self.theta) * cos(self.psi) * self.u + ( sin(self.phi) * sin(self.theta) * cos(self.psi) - cos(self.phi) * sin(self.psi)) * self.v + (cos(self.phi) * sin(self.theta) * cos(self.psi) + sin(self.phi)*sin(self.psi))*self.w ;
                    dot_pe      = cos(self.theta) * sin(self.psi) * self.u + ( sin(self.phi) * sin(self.theta) * sin(self.psi) + cos(self.phi) * cos(self.psi)) * self.v + (cos(self.phi) * sin(self.theta) * sin(self.psi) - sin(self.phi)*cos(self.psi))*self.w ;
                    dot_pd      = - sin(self.theta) * self.u + self.v * sin(self.phi) * cos(self.theta) + self.w * cos(self.phi) * cos(self.theta);
                    dot_u       = self.r * self.v - self.q * self.w - self.gravity * sin(self.theta) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * (CXalpha + CXqalpha * (self.c * self.q)/(2 * self.Va) + CXdeltae * u1) + self.rho * self.S_prop * self.C_prop / (2 * self.m) * ((self.k_motor * u2)^2 - self.Va^2);
                    dot_v       = self.p * self.w - self.r * self.u + self.gravity * cos(self.theta) * sin(self.phi) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * (self.C_Y0 + self.C_Ybeta * self.beta + self.C_Yp * self.b * self.r / (2* self.Va) + self.C_Ydeltaa * u3 + self.C_Ydeltar * u4) ;
                    dot_w       = self.q * self.u - self.p * self.v + self.gravity * cos(self.theta) * cos(self.phi) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * ( CZalpha + CZqalpha *(self.c * self.q)/(2*self.Va) + CZdeltae *u1) ;
                    dot_phi     = self.p + self.q * sin(self.phi) * tan(self.theta) + self.r * cos(self.phi) * tan(self.theta);
                    dot_theta  = self.q * cos (self.phi) - self.r * sin(self.phi);
                    dot_psi     = self.q * sin(self.phi) * sec(self.theta) + self.r * cos(self.phi) * sec(self.theta);
                    dot_p        = self.GAMMA1 * self.p * self.q - self.GAMMA2 * self.q * self.r + 1/2 * self.rho * (self.Va)^2 * self.S * self.b * (self.C_p0 + self.C_pbeta * self.beta + self.C_pp * self.b * self.p / (2 * self.Va) + self.C_pr * self.b * self.r /(2*self.Va) + self.C_pdeltaa * u3 + self.C_pdeltar * u4)  ;
                    dot_q        = self.GAMMA5 * self.p * self.r - self.GAMMA6 * (self.p^2 - self.r^2) + self.rho * (self.Va)^2 * self.S * self.c / (2 * self.J_y) * (self.C_m0 + self.C_malpha * self.alpha + self.C_mdeltae * u1);
                    dot_r         = self.GAMMA7 * self.p * self.q - self.GAMMA1 * self.q * self.r + 1/2 * self.rho * (self.Va)^2 * self.S * self.b * ( self.C_r0 + self.C_rbeta * self.beta + self.C_rp* self.b * self.p /(2 * self.Va ) + self.C_rr * self.b * self.r / (2 * self.Va) + self.C_rdeltaa * u3 + self.C_rdeltar * u4 ) ;
                   
                % update the state
                    self.pn         = self.pn + dot_pn * u5;
                    self.pe         = self.pe + dot_pe * u5;
                    self.pd         = self.pd + dot_pd * u5;
                    self.u           = self.u + dot_u * u5;
                    self.v           = self.v + dot_v * u5;
                    self.w          = self.w + dot_w * u5;
                    self.phi        = self.phi + dot_phi * u5;
                    self.theta     = self.theta + dot_theta * u5;
                    self.psi         = self.psi + dot_psi * u5;
                    self.p           = self.p + dot_p * u5;
                    self.q           = self.q + dot_q * u5;
                    self.r           = self.r + dot_r * u5;

                    self.Va         = sqrt(self.u^2 + self.v^2 + self.w^2);
                    self.alpha     = atan(self.w/self.u);
                    self.beta       = asin(self.v/self.Va);
        end 

        function f = find_f(self, u1, u2, u3, u4)

                 % dynamic coef
                   CLalpha     = self.C_L0 + self.C_Lalpha * self.alpha;
                    CDalpha     = self.C_D0 + self.C_Dalpha *self.alpha;
                    CXalpha     = - CDalpha * cos(self.alpha) + CLalpha * sin(self.alpha);
                    CXqalpha   = - self.C_Dq * cos(self.alpha) + self.C_Lq * sin(self.alpha);
                    CXdeltae    = - self.C_Ddeltae * cos(self.alpha) + self.C_Ldeltae * sin(self.alpha);
                    CZalpha     = - CDalpha * sin(self.alpha) - CLalpha * cos(self.alpha);
                    CZqalpha   = - self.C_Dq * sin(self.alpha) - self.C_Lq * cos(self.alpha);
                    CZdeltae    = - self.C_Ddeltae * sin(self.alpha) - self.C_Ldeltae * cos(self.alpha);

                % deriaviate of state
                    dot_pd      = - sin(self.theta) * self.u + self.v * sin(self.phi) * cos(self.theta) + self.w * cos(self.phi) * cos(self.theta);
                    dot_u       = self.r * self.v - self.q * self.w - self.gravity * sin(self.theta) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * (CXalpha + CXqalpha * (self.c * self.q)/(2 * self.Va) + CXdeltae * u1) + self.rho * self.S_prop * self.C_prop / (2 * self.m) * ((self.k_motor * u2)^2 - self.Va^2);
                    dot_v       = self.p * self.w - self.r * self.u + self.gravity * cos(self.theta) * sin(self.phi) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * (self.C_Y0 + self.C_Ybeta * self.beta + self.C_Yp * self.b * self.p / (2* self.Va) + self.C_Yr  * self.b * self.r /(2* self.Va) + self.C_Ydeltaa * u3 + self.C_Ydeltar * u4) ;
                    dot_w       = self.q * self.u - self.p * self.v + self.gravity * cos(self.theta) * cos(self.phi) + self.rho * (self.Va)^2 * self.S /(2 * self.m) * ( CZalpha + CZqalpha *(self.c * self.q)/(2*self.Va) + CZdeltae *u1) ;
                    dot_phi     = self.p + self.q * sin(self.phi) * tan(self.theta) + self.r * cos(self.phi) * tan(self.theta);
                    dot_theta  = self.q * cos(self.phi) - self.r * sin(self.phi);
                    dot_psi     = self.q * sin(self.phi) * sec(self.theta) + self.r * cos(self.phi) * sec(self.theta);
                    dot_p        = self.GAMMA1 * self.p * self.q - self.GAMMA2 * self.q * self.r + 1/2 * self.rho * (self.Va)^2 * self.S * self.b * (self.C_p0 + self.C_pbeta * self.beta + self.C_pp * self.b * self.p / (2 * self.Va) + self.C_pr * self.b * self.r /(2*self.Va) + self.C_pdeltaa * u3 + self.C_pdeltar * u4)  ;
                    dot_q        = self.GAMMA5 * self.p * self.r - self.GAMMA6 * (self.p^2 - self.r^2) + self.rho * (self.Va)^2 * self.S * self.c / (2 * self.J_y) * (self.C_m0 + self.C_malpha * self.alpha + self.C_mdeltae * u1);
                    dot_r         = self.GAMMA7 * self.p * self.q - self.GAMMA1 * self.q * self.r + 1/2 * self.rho * (self.Va)^2 * self.S * self.b * ( self.C_r0 + self.C_rbeta * self.beta + self.C_rp* self.b * self.p /(2 * self.Va ) + self.C_rr * self.b * self.r / (2 * self.Va) + self.C_rdeltaa * u3 + self.C_rdeltar * u4 ) ;
                    
                    f = [dot_pd; dot_u; dot_v; dot_w; dot_phi; dot_theta; dot_psi; dot_p; dot_q; dot_r];
        end
    end


end
