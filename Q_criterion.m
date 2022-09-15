function Q = Q_criterion(X,Y,z,U,V,W)

    [DUDX,DUDY,DUDZ] = gradHO_3D(X,Y,z,U);
    [DVDX,DVDY,DVDZ] = gradHO_3D(X,Y,z,V);
    [DWDX,DWDY,DWDZ] = gradHO_3D(X,Y,z,W);

    Q = DUDX.*DVDY + DUDX.*DWDZ + DVDY.*DWDZ ...
      - DUDY.*DVDX - DUDZ.*DWDX - DVDZ.*DWDY;

end