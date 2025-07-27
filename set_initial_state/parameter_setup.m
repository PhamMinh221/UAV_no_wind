function [pn, pe, pd, pn_center, pe_center, direction, Va] = parameter_setup()
    
    % set pn
    while true
        pn = input('Initial value of pn: ');
        if isnumeric(pn) && isscalar(pn)
            break;
        end
    end

    % set pe
    while true
        pe = input('Initial value of pe: ');
        if isnumeric(pe) && isscalar(pe)
            break;
        end
    end

    % set pd
    while true
        pd = input('Initial value of pd: ');
        if isnumeric(pd) && isscalar(pd)
            break;
        end
    end

    % set pn_center
    while true
        pn_center = input('Initial value of pn of the center point: ');
        if isnumeric(pn_center) && isscalar(pn_center)
            break;
        end
    end

    % set pe_center
    while true
        pe_center = input('Initial value of pe of the center point: ');
        if isnumeric(pe_center) && isscalar(pe_center)
            break;
        end
    end

    while true
        direction = input('Choose the direction of motion: 1 for clockwise or -1 for anticlockwise: ');
        if isnumeric(pe_center) && isscalar(pe_center) && (direction == 1 || direction == -1)
            break;
        end
    end

    while true
        Va = input('Desired velocity: ');
        if isnumeric(Va) && isscalar(Va)
            break;
        end
    end
    
end
