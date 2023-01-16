function vjk =  compute_induced_velocity(colp,hs1,hs2)
    rj = colp-hs1;
    rk = colp-hs2;
    l = hs2-hs1;
    term1 = 1/(4*pi);
    term2 = (cross(rj,rk))/(norm(cross(rj,rk)))^2;
    term3 = dot(l,rj)/norm(rj);
    term4 = dot(l,rk)/norm(rk);
    vjk = term1*term2*(term3-term4);
end