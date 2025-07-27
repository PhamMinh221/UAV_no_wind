classdef animation_craft < handle
    %
    %    Create spacecraft animation
    %
    %--------------------------------
    properties
        body_handle
    	Vertices
    	Faces
    	facecolors
        plot_initialized
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = animation_craft
            self.body_handle = [];
            [self.Vertices, self.Faces, self.facecolors] = self.define_spacecraft();
            self.plot_initialized = 0;           
        end
        %---------------------------
        function self=update(self, state)
            if self.plot_initialized==0
                figure(1); clf;
                self=self.drawBody(state.pn, state.pe, state.pd,...
                                   state.phi, state.theta, state.psi);
                title('Spacecraft')
                xlabel('East')
                ylabel('North')
                zlabel('High')
                %view(32,47)  % set the vieew angle for figure
                view(0,90)
                axis([-10,10,-10,10,-10,10]);
                hold on
                grid on
                self.plot_initialized = 1;
            else
                self=self.drawBody(state.pn, state.pe, state.pd,... 
                                   state.phi, state.theta, state.psi);

            end
        end
        %---------------------------
        function self = drawBody(self, pn, pe, pd, phi, theta, psi)
            Vertices = self.rotate(self.Vertices, phi, theta, psi);   % rotate rigid body  
            Vertices = self.translate(Vertices, pn, pe, pd);     % translate after rotation
            % transform vertices from NED to ENU (for matlab rendering)
            R = [...
                0, 1, 0;...
                1, 0, 0;...
                0, 0, 1;...
                ];
            Vertices = R*Vertices;
            if isempty(self.body_handle)
                self.body_handle = patch('Vertices', Vertices', 'Faces', self.Faces,...
                                             'FaceVertexCData',self.facecolors,...
                                             'FaceColor','flat');
            else
                set(self.body_handle,'Vertices',Vertices','Faces',self.Faces);
                axis([-50 50 -50 50 -50 50])
                drawnow
            end
        end 
        %---------------------------
        function pts=rotate(self, pts, phi, theta, psi)
            % define rotation matrix (right handed)
            R_roll = [...
                        1, 0, 0;...
                        0, cos(phi), sin(phi);...
                        0, -sin(phi), cos(phi)];
            R_pitch = [...
                        cos(theta), 0, -sin(theta);...
                        0, 1, 0;...
                        sin(theta), 0, cos(theta)];
            R_yaw = [...
                        cos(psi), sin(psi), 0;...
                        -sin(psi), cos(psi), 0;...
                        0, 0, 1];
            R = R_roll*R_pitch*R_yaw;   % inertial to body
            R = R';  % body to inertial
            % rotate vertices
            pts = R*pts;
        end
        %---------------------------
        % translate vertices by pn, pe, pd
        function pts = translate(self, pts, pn, pe, pd)
            pts = pts + repmat([pn;pe;-pd],1,size(pts,2));
        end
        %---------------------------
        function [V, F, colors] = define_spacecraft(self)
            
            fuse_h                  =   1.5;
            fuse_w                 =    1.5;
            fuse_l1                 =   2;
            fuse_l2                 =   1;  
            fuse_l3                 =   5;
            wing_l                  =   1;
            wing_w                 =  6;
            tailwing_l              =   1;
            tailwing_w            = 4.5;
            tail_h                   =  1.5;
            
            %Define_vertices
            point1              =   [fuse_l1, 0, 0];
            point2              =   [fuse_l2, fuse_l2, -fuse_h/2 ];
            point3              =   [fuse_l2, -fuse_l2, -fuse_h/2 ];
            point4              =   [fuse_l2, -fuse_l2, fuse_h/2 ];
            point5              =   [fuse_l2, fuse_l2, fuse_h/2 ];
            point6              =   [-fuse_l3, 0, 0 ];
            point7              =   [0, wing_w/2, 0 ];
            point8              =   [-wing_l, wing_w/2, 0];
            point9              =   [-wing_l, -wing_w/2, 0 ];
            point10            =   [0, -wing_w/2, 0];
            point11            =   [-fuse_l3 + tailwing_l, tailwing_w/2, 0  ];
            point12            =   [-fuse_l3, tailwing_w/2, 0 ];
            point13            =   [-fuse_l3, -tailwing_w/2, 0  ];
            point14            =   [-fuse_l3 + tailwing_l, -tailwing_w/2, 0  ];
            point15            =   [-fuse_l3 + tailwing_l, 0, 0 ];
            point16            =   [-fuse_l3, 0, -tail_h  ];
           
            
            V = [point1; point2; point3; point4; point5; point6; point7; point8; point9; point10; point11; point12; point13; point14; point15; point16]';

            % define faces as a list of vertices numbered above
            F = [1, 2, 3;
                   2, 3, 4;
                   1, 4, 5;
                   1, 2, 5;
                   2, 5, 6;
                   3, 4, 6;
                   2, 3, 6;
                   4, 5, 6;
                   7, 8, 9;
                   8, 9, 10;
                   11, 12, 13;
                   13, 12, 14;
                   15, 16, 6;];

            % define colors for each face    
            %myred = [1, 0, 0];
            mygreen = [0, 1, 0];
            myblue = [0, 0, 1];
            myyellow = [1, 1, 0];
            mycyan = [0, 1, 1];

            colors = [myblue; myblue; myblue; myblue; myyellow; myyellow; myyellow; myyellow;mygreen; mygreen; mygreen; mygreen; mycyan;];

        end
    end
end