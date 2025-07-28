
function psi = set_psi(pn_desired, pe_desired, pn_center, pe_center, direction)
    delta_y = pn_desired - pn_center;
    delta_x = pe_desired - pe_center;
    psi = pi/2 - atan(delta_y/delta_x) + pi * sign(pe_center - pe_desired)*(1+sign(pe_center - pe_desired))/2 + pi/2 * direction ;
    while (psi >= 2*pi)
        psi = psi -2*pi;
    end

end
