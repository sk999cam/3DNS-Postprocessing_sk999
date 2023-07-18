function s = tecplot2structured(blk, gas, bcs, tecfile, zonlbl)
    % tecfile = '/Data/ojj23/dns-work/channel-flows/RANS/M1-2_Re2k2-bsl/tecplot.dat';
    [zone, VARlist] = tec2mat(tecfile,'debug');

    gas.rgas = gas.cp*(gas.gam-1)/gas.gam;
    s = RANSSlice(blk, gas, bcs);
    s.blk = blk;
    s.gas = gas;
    gam  = gas.gam;

    nz = strcmpi(zonlbl, strtrim({zone.title}));
    data = zone(nz).data;
    xv = data(strcmp('X',VARlist)).data';
    yv = data(strcmp('Y',VARlist)).data';
    rv = data(strcmp('RHO',VARlist)).data';
    uv = data(strcmp('U',VARlist)).data';
    vv = data(strcmp('V',VARlist)).data';
    wv = data(strcmp('W',VARlist)).data';
    pv = data(strcmp('P',VARlist)).data';
    kv = data(strcmp('K',VARlist)).data';
    ov = data(strcmp('O',VARlist)).data';

    ri = scatteredInterpolant(xv, yv, rv);
    ui = scatteredInterpolant(xv, yv, uv);
    vi = scatteredInterpolant(xv, yv, vv);
    wi = scatteredInterpolant(xv, yv, wv);
    pi = scatteredInterpolant(xv, yv, pv);
    ki = scatteredInterpolant(xv, yv, kv);
    oi = scatteredInterpolant(xv, yv, ov);

    

    for ib = 1:length(blk.x)
        fprintf('Interpolating block %d/%d\n',[ib, length(blk.x)]);
        disp('Interpolating ro')
        s.ro{ib} = ri(blk.x{ib}, blk.y{ib});
        disp('Interpolating u')
        s.u{ib} = ui(blk.x{ib}, blk.y{ib});
        disp('Interpolating v')
        s.v{ib} = vi(blk.x{ib}, blk.y{ib});
        disp('Interpolating w')
        s.w{ib} = wi(blk.x{ib}, blk.y{ib});
        disp('Interpolating p')
        pnow = pi(blk.x{ib}, blk.y{ib});
        s.Et{ib} = pnow/(gam-1) + 0.5*s.ro{ib}.*(s.u{ib}.^2 + s.v{ib}.^2 + s.w{ib}.^2);
        disp('Interpolating k')
        s.k{ib} = ki(blk.x{ib}, blk.y{ib});
        disp('Interpolating omega')
        s.omega{ib} = oi(blk.x{ib}, blk.y{ib});

    end
    s.getBCs;
end