function blparams = integralparams(x,y,qbar,mu_ref,gam,rgas,vref,nb,ns,xLE,xTE,blparams)

[ni,nj]=size(x);
cax = xTE-xLE;

ISS = blparams.ISS;
IPS = blparams.IPS;
XX = blparams.XX;
YY = blparams.YY;
SS = blparams.SS;

x0 = blparams.x0;
y0 = blparams.y0;

d0 = blparams.d0;

d1 = blparams.d1;
d2 = blparams.d2;
d3 = blparams.d3;

d1i = blparams.d1i;
d2i = blparams.d2i;
d3i = blparams.d3i;

d4 = blparams.d4;
d5 = blparams.d5;
%blparams.d6 = d6;
d7 = blparams.d7;
%blparams.d8 = d8;
%blparams.d9 = d9;
d10 = blparams.d10;
d11 = blparams.d11;
d12 = blparams.d12;
d13 = blparams.d13;
%blparams.d14 = d14;
d15 = blparams.d15;
d16 = blparams.d16;
d17 = blparams.d17;
ue = blparams.ue;
d3u3 = blparams.d3u3;
skinfriction = blparams.skinfriction;
retheta = blparams.retheta;
tke_edge = blparams.tke_edge;
pwall = blparams.pwall;

ro = qbar.ro;
ru = qbar.ru;
rv = qbar.rv;
rw = qbar.rw;
Et = qbar.Et;

ruu = qbar.ruu;
rvv = qbar.rvv;
rww = qbar.rww;

ruv = qbar.ruv;
ruw = qbar.ruw;
rvw = qbar.rvw;

%rus = qbar.rus;
%diss = qbar.diss;

u = ru./ro;
v = rv./ro;
w = rw./ro;
p = (gam-1)*(Et - 0.5*(ruu + rvv + rww));
T = p./(ro*rgas);

% s = cp*log(T) - rgas*log(p);
% 
%  s = rus./ru;
%  ii = isnan(s);
%  s(ii)=0;

vel = sqrt(u.*u + v.*v + w.*w);
mach = vel./sqrt(gam*rgas*T);
To = T.*(1.0 + (gam-1)*0.5*mach.*mach);
po = p.*((To./T).^(gam/(gam-1)));

mu = (mu_ref)*( ( 273.0 + 110.4 )./( T + 110.4 ) ).*((T/273).^1.5) ;
%cond = mu*cp/0.72;

%dissT = diss./T;

[dudx,dudy] = gradHO(x,y,u);
[dvdx,dvdy] = gradHO(x,y,v);
[dwdx,dwdy] = gradHO(x,y,w);

s11 = dudx;
s22 = dvdy;
s33 = -s11-s22;

s12 = 0.5*(dudy + dvdx);
s23 = 0.5*dwdy;
s13 = 0.5*dwdx;

dwdz = s33;

diss_av = ((mu.*(2*(s11.*s11 + s22.*s22 + s33.*s33) + 4*s23.*s23 + 4*s13.*s13 + 4*s12.*s12 ) ...
            - (2/3)*mu.*(s11 + s22 + s33).*(s11 + s22 + s33) ));
diss_avT = diss_av./T;

ruvdash =(ruv - ru.*v);
ruwdash =(ruw - ru.*w);
rvwdash =(rvw - rv.*w);
ruudash =(ruu - ru.*u);
rvvdash =(rvv - rv.*v);
rwwdash =(rww - rw.*w);
 
turb =  -( ruvdash.*(dudy + dvdx) + rvwdash.*dwdy + ruwdash.*dwdx ...
          + ruudash.*dudx + rvvdash.*dvdy + rwwdash.*dwdz);

turbT = turb./T;
eke = (ruudash + rvvdash + rwwdash)*0.5./ro;

%misen = sqrt((((p./poin).^(-(gam-1)/gam))-1)*2/(gam-1));
%visen = misen.*sqrt(gam*rgas*T);
%ke_isen = 0.5*visen.*visen;
%ke_defect = ke_isen - 0.5*(ruu + rvv + rww)./ro + eke;

[dtkedx,dtkedy] = gradHO(x,y,eke);
%[dkedefdx,dkedefdy] = gradHO(x,y,ke_defect);

div_tke = u.*dtkedx + v.*dtkedy;
%div_kedef = u.*dkedefdx + v.*dkedefdy;

[dpdx,dpdy] = gradHO(x,y,p);
ke_isen_flux = -(u.*dpdx + v.*dpdy)./ro;

kebar = 0.5*(ruu./ro + rvv./ro + rww./ro);

[dkedx,dkedy] = gradHO(x,y,0.5*(u.*u + v.*v + w.*w));

ke_flux = (u.*dkedx + v.*dkedy);%-div_tke;

diss_act = (turb./ro - div_tke);


istart = 1;
iend = ni;
istep = 1;

if nb==7
istart = ni;
iend = 1;
istep = -1;    
end


ii = (istart + sign(istep)):istep:(iend-sign(istep));

dxt(ii,1) = (x(ii+sign(istep),nj,1)- x(ii-sign(istep),nj,1));
dyt(ii,1) = (y(ii+sign(istep),nj,1)- y(ii-sign(istep),nj,1));

dxt(iend,1) = ((x(iend,nj,1)- x(iend-sign(istep),nj,1)))*2.0;
dyt(iend,1) = ((y(iend,nj,1)- y(iend-sign(istep),nj,1)))*2.0;

dxt(istart,1) = ((x(istart+sign(istep),nj,1)- x(istart,nj,1)))*2.0;
dyt(istart,1) = ((y(istart+sign(istep),nj,1)- y(istart,nj,1)))*2.0;

dxt = 0.5*dxt;
dyt = 0.5*dyt;

wall_dist = sqrt( (x(:,nj-1,1)- x(:,nj,1)).^2 + (y(:,nj-1,1) - y(:,nj,1)).^2 );
tang_dist = sqrt(dxt.*dxt + dyt.*dyt);
uwall = (u(:,nj-1).*dxt + v(:,nj-1).*dyt)./sqrt(dxt.*dxt + dyt.*dyt);

tauw = mu(:,nj).*uwall./wall_dist;
utau = sqrt(abs(tauw)./ro(:,nj));

% yplus = utau.*ro(:,nj).*wall_dist./mu(:,nj);
% xplus = yplus.*tang_dist./wall_dist;
% zplus = (span(2)-span(1))*yplus./wall_dist;

minf = sqrt((((p(:,nj)./max(po,[],2)).^(-(gam-1)/gam))-1)*2/(gam-1));
%Tinf = To(:,1).*((p(:,nj)./max(po,[],2)).^((gam-1)/gam));
Tinf = To(:,nj).*((p(:,nj)./max(po,[],2)).^((gam-1)/gam));

uinf = minf.*sqrt(gam*rgas*Tinf);
vinf = vref;
uisen = uinf;

uinf2 = max(sqrt(u.*u + v.*v + w.*w),[],2);
%uinf2 = max(sqrt(u.*u + v.*v),[],2);

ie = uinf > uinf2;
uinf(ie) = uinf2(ie);

roinf = p(:,nj)./(rgas*Tinf);

ruinf = roinf.*uinf;
cf = tauw./(0.5*ruinf.*uinf);

for i=1:ni

ns = ns + 1;  

% guess b.layer height
if(ns>1)
hprev = d1(ns-1)/d2(ns-1);
del = d2(ns-1).*(3.15 + 1.72./(hprev-1)) + d1(ns-1);    
else
del = 1e-4;
end

ydist = y(i,nj-3)-y(i,nj);
xdist = x(i,nj-3)-x(i,nj);
% 
IPS(ns) = ydist < 0 & (x(i,nj) < (xTE-cax*0.01)) ;%| xdist > 0;
ISS(ns) = ydist > 0 & (x(i,nj) < (xTE-cax*0.01));%| xdist < 0;

volbl(ns) = 0;

x0(ns) = 0;
y0(ns) = 0;
d0(ns) = 0;

d1(ns) = 0;
d2(ns) = 0;
d3(ns) = 0;
d1i(ns) = 0;
d2i(ns) = 0;
d3i(ns) = 0;


d4(ns) = 0;
d5(ns) = 0;
d6(ns) = 0;
d7(ns) = 0;
d8(ns) = 0;
d9(ns) = 0;
d10(ns) = 0;
d11(ns) = 0;
d12(ns) = 0;
d13(ns) = 0;
d14(ns) = 0;
d15(ns) = 0;
d16(ns) = 0;
d17(ns) = 0;
d18(ns) = 0;
d19(ns) = 0;
d20(ns) = 0;
d21(ns) = 0;
d22(ns) = 0;
d23(ns) = 0;
d24(ns) = 0;

iedge = 0;

ynow = 0;
j0(ns) = nj;


yprof = y(i,:);
xprof = x(i,:);
uprof = u(i,:);
vprof = v(i,:);
wprof = w(i,:);

roprof = ro(i,:);
Ekprof = eke(i,:);
%dissprof = diss(i,:)./ro(i,:);
%dissprof = diss_act(i,:);
%dissTprof = dissT(i,:);
turbprof = turb(i,:)./ro(i,:);
keprof = kebar(i,:);
divtkeprof = div_tke(i,:);
%divkedefprof = div_kedef(i,:);
dissavprof = diss_av(i,:)./ro(i,:);
uvprof = ruvdash(i,:)./ro(i,:);
%ke_isen_flux_prof = ke_isen_flux(i,:);
ke_flux_prof = ke_flux(i,:);
velprof = vel(i,:);

nprof = 10000;

xd = xprof-xprof(nj);
yd = yprof-yprof(nj);
rprof = sqrt(xd.*xd + yd.*yd);

ri = linspace(min(rprof),max(rprof),nprof);
yi = interp1(rprof,yprof,ri,'pchip');
xi = interp1(rprof,xprof,ri,'pchip');
ui = interp1(rprof,uprof,ri,'pchip');
vi = interp1(rprof,vprof,ri,'pchip');
wi = interp1(rprof,wprof,ri,'pchip');

roi = interp1(rprof,roprof,ri,'pchip');
Eki = interp1(rprof,Ekprof,ri,'pchip');
%dissi = interp1(rprof,dissprof,ri,'pchip');
%dissTi = interp1(rprof,dissTprof,ri,'pchip');
turbi = interp1(rprof,turbprof,ri,'pchip');
kebari = interp1(rprof,keprof,ri,'pchip');
divtkei = interp1(rprof,divtkeprof,ri,'pchip');
%divkedefi = interp1(rprof,divkedefprof,ri,'pchip');
dissavi = interp1(rprof,dissavprof,ri,'pchip');
uvi = interp1(rprof,uvprof,ri,'pchip');
%ke_isen_fluxi = interp1(rprof,ke_isen_flux_prof,ri,'pchip');
ke_fluxi = interp1(rprof,ke_flux_prof,ri,'pchip');

vmagi = sqrt(ui.*ui + vi.*vi + wi.*wi);

iedge = 0;

for n=2:nprof
% 

if(ns<100)
 if(iedge==0 && vmagi(n) >   uinf(i)*0.99999 )
   iedge = 1;
   tke_edge(ns) = Eki(n);
   [~,edgej(i)] = min(abs(rprof-ri(n)));
 end
else
 if(iedge==0 && d0(ns) > del )
  iedge = 1;
  tke_edge(ns) = Eki(n);
  [~,edgej(i)] = min(abs(rprof-ri(n)));
 end
end

dx = (xi(n) - xi(n-1));
dy = (yi(n) - yi(n-1));

dely = abs(-dx*dyt(i) + dy*dxt(i))/sqrt(dxt(i)*dxt(i) + dyt(i)*dyt(i)); 
dveldy = (vmagi(n)-vmagi(n-1))/dely;

uav = (ui(n) + ui(n-1))*0.5;
vav = (vi(n) + vi(n-1))*0.5;
wav = (wi(n) + wi(n-1))*0.5;

Ekav = (Eki(n) + Eki(n-1))*0.5;
%dissav = (dissi(n) + dissi(n-1))*0.5;
%dissTav = (dissTi(n) + dissTi(n-1))*0.5;
turbav = (turbi(n) + turbi(n-1))*0.5;
keav = (kebari(n) + kebari(n-1))*0.5 - Ekav;
divtkeav = (divtkei(n) + divtkei(n-1))*0.5;
%divkedefav = (divkedefi(n) + divkedefi(n-1))*0.5;
dissavav = (dissavi(n) + dissavi(n-1))*0.5;
uvav = (uvi(n) + uvi(n-1))*0.5;
%ke_isen_fluxav = (ke_isen_fluxi(n) + ke_isen_fluxi(n-1))*0.5;
ke_fluxav = (ke_fluxi(n) + ke_fluxi(n-1))*0.5;

ronow = (roi(n) + roi(n-1))*0.5;
vnow = (uav*dxt(i) + vav*dyt(i))/sqrt(dxt(i)*dxt(i) + dyt(i)*dyt(i)); 

if(vnow>uinf(i))
    vnow = uinf(i);
end

ynow = ynow+dely;
if(iedge~=1)
x0(ns) = (xi(n) + xi(n-1))*0.5;
y0(ns) = (yi(n) + yi(n-1))*0.5;
end

d0(ns) = d0(ns) + dely*(1-iedge);
d1(ns) = d1(ns) + ( 1 - ronow.*vnow/ruinf(i) )*dely*(1-iedge);
%d1(ns) = d1(ns) + ( 1 - vnow/uinf(i) )*dely*(1-iedge);
d2(ns) = d2(ns) + (ronow.*vnow/ruinf(i))*( 1 - vnow/uinf(i) )*dely*(1-iedge);   
d3(ns) = d3(ns) + (ronow.*vnow/ruinf(i))*( 1 - (vnow/uinf(i))*(vnow/uinf(i)) )*dely*(1-iedge);   
d4(ns) = d4(ns) + vnow*Ekav*dely*(1-iedge);   

d1i(ns) = d1i(ns) + ( 1 - vnow/uinf(i) )*dely*(1-iedge);
d2i(ns) = d2i(ns) + (vnow/uinf(i))*( 1 - vnow/uinf(i) )*dely*(1-iedge);   
d3i(ns) = d3i(ns) + (vnow/uinf(i))*( 1 - (vnow/uinf(i))*(vnow/uinf(i)) )*dely*(1-iedge);   


%d5(ns) = d5(ns) + dissav*dely*(1-iedge);   
%d6(ns) = d6(ns) + efluxav*dely*(1-iedge);   
%d7(ns) = d7(ns) + dissTav*dely*(1-iedge);
%d8(ns) = d8(ns) + rdsdtav*dely*(1-iedge);
%d9(ns) = d9(ns) + thermav*dely*(1-iedge);
d10(ns) = d10(ns) + turbav*dely*(1-iedge);
d11(ns) = d11(ns) + divtkeav*dely*(1-iedge);
%d12(ns) = d12(ns) + divkedefav*dely*(1-iedge);
d13(ns) = d13(ns) + dissavav*dely*(1-iedge);
%d14(ns) = d14(ns) + diffav*dely*(1-iedge);
d15(ns) = d15(ns) + vnow*vnow*0.5*dely*(1-iedge);
%d16(ns) = d16(ns) + ke_isen_fluxav*dely*(1-iedge);
d17(ns) = d17(ns) + ke_fluxav*dely*(1-iedge);
% d18(ns) = d18(ns) + visav*dely*(1-iedge);
% d19(ns) = d19(ns) + turb2av*dely*(1-iedge);
% d20(ns) = d20(ns) + turb3av*dely*(1-iedge);
% d21(ns) = d21(ns) + norm2av*dely*(1-iedge);
% d22(ns) = d22(ns) + normav*dely*(1-iedge);
% d23(ns) = d23(ns) + shearav*dely*(1-iedge);
% d24(ns) = d24(ns) + efluxavav*dely*(1-iedge);
end

d3u3(ns) = d3(ns)*(uinf(i)^3);
retheta(ns) = d2(ns)*ruinf(i)/mu(i,nj);
u3(ns) = (uinf(i)^3);
ue(ns) = uinf(i);
skinfriction(ns) = cf(i);
if(i>1)
snow(ns) = snow(ns-1) + tang_dist(i);
else
snow(ns)=0;
end

XX(ns) = x(i,nj);
YY(ns) = y(i,nj);

pwall(ns) = p(i,nj);

end

blparams.x0 = x0;
blparams.y0 = y0;

blparams.ISS = ISS;
blparams.IPS = IPS;
blparams.XX = XX;
blparams.YY = YY;
blparams.SS = snow;

blparams.d0 = d0;

blparams.d1 = d1;
blparams.d2 = d2;
blparams.d3 = d3;

blparams.d1i = d1i;
blparams.d2i = d2i;
blparams.d3i = d3i;

blparams.d4 = d4;
blparams.d5 = d5;
%blparams.d6 = d6;
blparams.d7 = d7;
%blparams.d8 = d8;
%blparams.d9 = d9;
blparams.d10 = d10;
blparams.d11 = d11;
blparams.d12 = d12;
blparams.d13 = d13;
%blparams.d14 = d14;
blparams.d15 = d15;
blparams.d16 = d16;
blparams.d17 = d17;
blparams.ue = ue;
blparams.d3u3 = d3u3;
blparams.skinfriction = skinfriction;
blparams.retheta = retheta;
blparams.tke_edge = tke_edge;
blparams.edgej = edgej;
blparams.pwall = pwall;

return