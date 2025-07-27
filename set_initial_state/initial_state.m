

% State variables for MAV equations of motion
    MAV.pn          = 0;            % Inertial North position 
    MAV.pe          = -40;            % Inertial East position 
    MAV.pd          = 0;            % Inertial Down position 

    MAV.u            = 17;            % Body frame velocity along i
    MAV.v            = 0;                       % Body frame velocity along j
    MAV.w            = 0;            % Body frame velocity along k

    MAV.phi          = 0;                       % Roll angle
    MAV.theta       = 0;            % Pitch angle
    MAV.psi          = 0;                   % Yaw angle

    MAV.p            = 0;             % Roll rate along i in Fb
    MAV.q            = 0;             % Pitch rate along j in Fb
    MAV.r             = 0;              % Yaw rate along k in Fb

    MAV.Va           =17;               % velocity relative to the wind
 % Other angle
    MAV.alpha      = 0;              % attack angle
    MAV.beta        = 0;              % side slip angle
    MAV.chi          = 0;               % course angle
    MAV.gamma    = 0;               % flight path angle
    MAV.chic         = MAV.chi - MAV.psi;       %crab angle
    MAV.gammaa  = MAV.theta - MAV.alpha;        %air mass refferenced flight path angle
    
  % derivative 
    
