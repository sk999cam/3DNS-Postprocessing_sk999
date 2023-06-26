function mu_ref = sutlerland_mu_ref(mu, T, mu_cref, mu_tref)
    if nargin < 3
        mu_cref = 110.4;
        mu_tref = 273;
    end

    mu_ref = mu * (mu_tref/T)^1.5 * (T + mu_cref) / (mu_tref + mu_cref);
end