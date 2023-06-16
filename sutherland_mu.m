function mu = sutherland_mu(T, mu_ref, mu_cref, mu_tref)
    if nargin < 3
        mu_cref = 110.4;
        mu_tref = 273;
    end
    if nargin < 2
        mu_ref = 1.716e-5;
    end

    mu = mu_ref * (T/mu_tref)^1.5 * (mu_tref + mu_cref) / (T + mu_cref);
    
end