function StR = strain_rate_magnitude(x, y, u, v)
    [DUDX, DUDY] = gradHO(x, y, u);
    [DVDX, DVDY] = gradHO(x, y, v);
    StR = sqrt(2*DUDX.^2 + 2*DVDY.^2 + (DUDY + DVDX).^2);
end