clear;clc;

% addpath(genpath('C:\Users\minh2\OneDrive\Desktop\UAV\Pham_Hong_Minh_VDT_Baitap_UAV'));
rootPath = fileparts(mfilename('fullpath'));
addpath(genpath(rootPath));

MAV = UAV;
run("parameter_aerosonde.m");
run("initial_state");
time_sampling = 0.01;


%User set up
[MAV.pn, MAV.pe, MAv.pd, pn_center, pe_center, direction, Va_command ] =  parameter_setup();
Radius_of_trajectory= sqrt((MAV.pn - pn_center)^2 + (MAV.pe - pe_center)^2);
MAV.psi = set_psi(MAV.pn, MAV.pe, pn_center,pe_center, direction);

[state_trim, input_trim] = solve_trim(MAV, Va_command, Radius_of_trajectory); 


MAV.u           = state_trim(1);
MAV.v           = state_trim(2);
MAV.w           = state_trim(3);
MAV.theta      = state_trim(4);
MAV.p           = state_trim(5);
MAV.q           = state_trim(6);
MAV.r            = state_trim(7);
MAV.alpha     = state_trim(8);
MAV.beta      = state_trim(9);
MAV.phi        = state_trim(10);
MAV.psi = set_psi(MAV.pn, MAV.pe, pn_center,pe_center, direction);


delta_e         = input_trim(1);
delta_t          = input_trim(2);
delta_a          = input_trim(3);
delta_r          = input_trim(4);

mav_simulation = animation_craft();

last_point = [MAV.pe, MAV.pn, MAV.pd];

for t = 1:2000

   MAV = MAV.update_state(delta_e, delta_t, delta_a, delta_r, time_sampling);
   tic
   mav_simulation.update(MAV);

   current_point = [MAV.pe, MAV.pn, MAV.pd];

   plot3([last_point(1), current_point(1)], ...
          [last_point(2), current_point(2)], ...
          -[last_point(3), current_point(3)], ...
          'r-', 'LineWidth', 1.2);
   last_point = current_point;
   hold on;

   drawnow limitrate
   runtime = toc;
   disp(runtime);
end



