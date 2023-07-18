function val = p_p0(M, gam)

if nargin<2
    gam=1.4;
end

val = (1 + 0.5*(gam-1)*M.^2).^(-gam/(gam-1));

end