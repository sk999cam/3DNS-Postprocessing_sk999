function [flow,geom,blparams,inlet,exit,wake,suction] = process_mean(dir)
% 
run([dir,'case_data'])

load([dir,'blockdims.txt']);

nk = blockdims(1,3);

span = linspace(0,span,nk);

%nk = 1;

cv = cp/gam;
rgas = cp-cv;

iface_tot = 0;
nn = 0;
ns = 0;
npoints = 0;
pstag = 0;

poin = ptot;

cax = xTE-xLE;
xwake = xTE + 0.5*cax;
xexit = xTE + 0.75*cax;%0.25*cax;
xinlet = xLE - 0.25*cax;
xblprof = xTE - 0.03*cax;

load([dir,'mean_time.txt']);
tot_time = sum(mean_time(nmean_start:nmean_end,3));

%nphase = 8;
nstats_mean = 17;
nstats = nstats_mean + 11*nphase; 

minx = 1e12;
maxx = -1e12;
miny = 1e12;
maxy = -1e12;

for nb=1:Nb

    
fid = fopen([dir,'grid_',num2str(nb),'.txt'],'r');

   ni = blockdims(nb,1);
   nj = blockdims(nb,2);
      
    x = zeros(ni,nj);
    y = zeros(ni,nj);
    
   
   for j=1:nj
   for i=1:ni
   [A]=fscanf(fid, '%f %f',2);
   x(i,j) = A(1);
   y(i,j) = A(2);
   end
   end
   
   fclose(fid);


xb{nb} = x;
yb{nb} = y;


geom{nb}.x = x;
geom{nb}.y = y;


end


nk=1;
for nb=1:Nb
 
   x = xb{nb};
   y = yb{nb};
   
for nmean=nmean_start:nmean_end
     
  
delt_now = mean_time(nmean,3);
     
   
   ni = blockdims(nb,1);
   nj = blockdims(nb,2);
 

   if(nmean==nmean_start)
   ro = zeros(ni,nj,nk);
   ru = zeros(ni,nj,nk);
   rv = zeros(ni,nj,nk);
   rw = zeros(ni,nj,nk);
   Et = zeros(ni,nj,nk);
   
   r2 =  zeros(ni,nj,nk);
   ruu =  zeros(ni,nj,nk);
   rvv =  zeros(ni,nj,nk);
   rww =  zeros(ni,nj,nk);
   
   ruv =  zeros(ni,nj,nk);
   ruw =  zeros(ni,nj,nk);
   rvw =  zeros(ni,nj,nk);
   
   p2 =  zeros(ni,nj,nk);
   
   rus =  zeros(ni,nj,nk);
   rvs =  zeros(ni,nj,nk);
   rws =  zeros(ni,nj,nk);
   
   diss = zeros(ni,nj,nk);
   snow = zeros(ni,nj,nk);
   pnow = zeros(ni,nj,nk);
   tnow = zeros(ni,nj,nk);
   ponow = zeros(ni,nj,nk);
   tonow = zeros(ni,nj,nk);
   end   
   
%

fid2 = fopen([dir,'mean2_',num2str(nb),'_',num2str(nmean)],'r');    
A = fread(fid2,inf,'float64');
A = reshape(A,nstats,length(A)/nstats);

fid3 = fopen([dir,'mnod2_',num2str(nb),'_',num2str(nmean)],'r');    
B = fread(fid3,inf,'int');
B = reshape(B,3,length(B)/3);


icount(1:ni*nj*3) = 0;

for n=1:size(A,2)
    
i = B(1,n);
j = B(2,n);
k = B(3,n);


nid = i + (j-1)*ni;

if(i<=ni && j<=nj && icount(nid)==0)
    
ro(i,j) = ro(i,j) + A(1,n)/tot_time;
ru(i,j) = ru(i,j) + A(2,n)/tot_time;
rv(i,j) = rv(i,j) + A(3,n)/tot_time;
rw(i,j) = rw(i,j) + A(4,n)/tot_time;
Et(i,j) = Et(i,j) + A(5,n)/tot_time;

r2(i,j) =   r2(i,j) + A(6,n)/tot_time;
ruu(i,j) = ruu(i,j) + A(7,n)/tot_time;
rvv(i,j) = rvv(i,j) + A(8,n)/tot_time;
rww(i,j) = rww(i,j) + A(9,n)/tot_time;

ruv(i,j) = ruv(i,j) + A(10,n)/tot_time;
ruw(i,j) = ruw(i,j) + A(11,n)/tot_time;
rvw(i,j) = rvw(i,j) + A(12,n)/tot_time;
p2(i,j) =   p2(i,j) + A(13,n)/tot_time;

rus(i,j) =   rus(i,j) + A(16,n)/tot_time;
diss(i,j) = diss(i,j) + A(17,n)/tot_time;

end

icount(nid) = 1;

end

fclose(fid2);
fclose(fid3);

end
    

u = ru./ro;
v = rv./ro;
w = rw./ro;
p = (gam-1)*(Et - 0.5*(ruu + rvv + rww));
T = p./(ro*rgas);

s = cp*log(T) - rgas*log(p);

 s = rus./ru;
 ii = isnan(s);
 s(ii)=0;

vel = sqrt(u.*u + v.*v + w.*w);
%vel = sqrt(ruu./ro + rvv./ro + rww./ro);

mach = vel./sqrt(gam*rgas*T);

To = T.*(1.0 + (gam-1)*0.5*mach.*mach);
po = p.*((To./T).^(gam/(gam-1)));

mu = (mu_ref)*( ( 273.0 + 110.4 )./( T + 110.4 ) ).*((T/273).^1.5) ;
cond = mu*cp/0.72;

dissT = diss./T;

% [drusdx,drusdy] = gradHO(x,y,rus);
% [drvsdx,drvsdy] = gradHO(x,y,rvs);
% [drudx,drudy] = gradHO(x,y,ru);
% [drvdx,drvdy] = gradHO(x,y,rv);
% [dTdx,dTdy] = gradHO(x,y,T);
% [d2Tx,~] = gradHO(x,y,dTdx);
% [~,d2Ty] = gradHO(x,y,dTdy);
% [drusdx,drusdy] = gradHO(x,y,rus);
% [drvsdx,drvsdy] = gradHO(x,y,rvs);
[drudx,drudy] = gradHO(x,y,ru);
[drvdx,drvdy] = gradHO(x,y,rv);
[dTdx,dTdy] = gradHO(x,y,T);
[d2Tx,~] = gradHO(x,y,dTdx);
[~,d2Ty] = gradHO(x,y,dTdy);

[dqx,~] = gradHO(x,y,-cond.*dTdx./T);
[~,dqy] = gradHO(x,y,-cond.*dTdy./T);

[dsavdx,dsavdy] = gradHO(x,y,s);


%ii = abs(ru) > 0.0;
%s(ii) = rus(ii)./ru(ii);

%rudsdx = drusdx - snow.*drudx;
%rvdsdy = drvsdy - snow.*drvdy;
%rudsdx = drusdx - s.*drudx;
%rvdsdy = drvsdy - s.*drvdy;
%rvgrads = rudsdx + rvdsdy;
rvgradsav = ru.*dsavdx + rv.*dsavdy;

rvgrads = rvgradsav;

%rvgrads = drusdx + drvsdy;

%rdsdt = ro.*(snow-sold)/sum(mean_time(nmean_start:nmean_end-1,3));
rdsdt = zeros(ni,nj);

therm = cond.*(d2Tx + d2Ty)./T;
%therm = -(dqx + dqy);


%ii = abs(ru) > 0.0;
%s(ii) = rus(ii)./ru(ii);

%s = (rus + rvs)./sqrt(ru.*ru + rv.*rv);


prms = sqrt(abs(p2 - p.*p));
urms = sqrt(abs(ruu./ro - u.*u));
vrms = sqrt(abs(rvv./ro - v.*v));
wrms = sqrt(abs(rww./ro - w.*w));

eke = 0.5*(urms.*urms + vrms.*vrms + wrms.*wrms);
tu = sqrt(eke*2/3)/vref;

% [drdx,drdy] = gradHO2(x,y,ro);
 [dudx,dudy] = gradHO(x,y,u);
 [dvdx,dvdy] = gradHO(x,y,v);
 [dwdx,dwdy] = gradHO(x,y,w);
 
 [d2udx,~] = gradHO(x,y,dudx);
 [d2vdx,~] = gradHO(x,y,dvdx);
 [d2wdx,~] = gradHO(x,y,dwdx);
 
 [~,d2udy] = gradHO(x,y,dudy);
 [~,d2vdy] = gradHO(x,y,dvdy);
 [~,d2wdy] = gradHO(x,y,dwdy);
 
% [drudx,~] = gradHO2(x,y,ru);
% [~,drvdy] = gradHO2(x,y,rv);
% [drusdx,~] = gradHO2(x,y,rus);
% [~,drvsdy] = gradHO2(x,y,rvs);
% [dsdx,dsdy] = gradHO2(x,y,s);

vis = mu.*( u.*(d2udx + d2udy) + v.*(d2vdx + d2vdy) + w.*(d2wdx + d2wdy) )./ro;

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

% 
% div_mass = drudx + drvdy;
% div_ent = drusdx + drvsdy;
% %div_s = ru.*dsdx + rv.*dsdy;% + s.*div_mass;
% div_s = u.*dsdx + v.*dsdy;

ruvdash =(ruv - ru.*v);
ruwdash =(ruw - ru.*w);
rvwdash =(rvw - rv.*w);
ruudash =(ruu - ru.*u);
rvvdash =(rvv - rv.*v);
rwwdash =(rww - rw.*w);


 [druvdx,druvdy] = gradHO(x,y,ruvdash);
 [druudx,~] = gradHO(x,y,ruudash);
 [~,drvvdy] = gradHO(x,y,rvvdash);
 
turb =  -( ruvdash.*(dudy + dvdx) + rvwdash.*dwdy + ruwdash.*dwdx ...
          + ruudash.*dudx + rvvdash.*dvdy + rwwdash.*dwdz);
norm =  -(ruudash.*dudx + rvvdash.*dvdy);% + rwwdash.*dwdz);
shear = -( ruvdash.*(dudy + dvdx) );%+ rvwdash.*dwdy + ruwdash.*dwdx );

[t1,~] = gradHO(x,y,u.*ruudash + v.*ruvdash + w.*ruwdash);
[~,t2] = gradHO(x,y,u.*ruvdash + v.*rvvdash + w.*rvwdash);
[t3,~] = gradHO(x,y,u.*ruudash);
[~,t4] = gradHO(x,y,v.*rvvdash);

turb2 = -(t1 + t2);
turb3 = -(v.*(drvvdy + druvdx) + u.*(druudx + druvdy)) ;
norm2 = -(t3 + t4);

turbT = turb./T;
eke = (ruudash + rvvdash + rwwdash)*0.5./ro;
misen = sqrt((((p./poin).^(-(gam-1)/gam))-1)*2/(gam-1));
visen = misen.*sqrt(gam*rgas*T);
ke_isen = 0.5*visen.*visen;

ke_defect = ke_isen - 0.5*(ruu + rvv + rww)./ro + eke;

[dtkedx,dtkedy] = gradHO(x,y,eke);
[dkedefdx,dkedefdy] = gradHO(x,y,ke_defect);

div_tke = u.*dtkedx + v.*dtkedy;
div_kedef = u.*dkedefdx + v.*dkedefdy;

diffx = -(mu./ro).*(dtkedx + (druudx./ro) + (druvdy./ro));
diffy = -(mu./ro).*(dtkedy + (druvdx./ro) + (drvvdy./ro));

[ddiffxdx,~] = gradHO(x,y,diffx);
[ddiffydy,~] = gradHO(x,y,diffy);

diffusion = ddiffxdx + ddiffydy;

thet = repmat(atan(v(:,nj-10)./u(:,nj-10)),[1 nj]);
ct2 = cos(thet).^2;
st2 = sin(thet).^2;
cs = 2.0*cos(thet).*sin(thet);

ruudashres = ruudash.*ct2 + rvvdash.*st2 + ruvdash.*cs;
rvvdashres = ruudash.*st2 + rvvdash.*ct2 - ruvdash.*cs;

pturb = ro.*eke.*2/3;
[dpdx,dpdy] = gradHO(x,y,p);
ke_isen_flux = -(u.*dpdx + v.*dpdy)./ro;

kebar = 0.5*(ruu./ro + rvv./ro + rww./ro);

[dkedx,dkedy] = gradHO(x,y,0.5*(u.*u + v.*v + w.*w));
[d2kedx,~] = gradHO(x,y,dkedx);
[~,d2kedy] = gradHO(x,y,dkedy);

ke_flux = (u.*dkedx + v.*dkedy);%-div_tke;

diss_act = (turb./ro - div_tke);
diffK = mu.*(d2kedx + d2kedy);

l_kolmog = (((mu./ro).^3.0)./abs(diss_act-diss_av)).^0.25;
% l_taylor = sqrt( 10*(mu./ro).*eke./abs(diss_tke./ro) );
% l_outer = (eke.*1.5)./abs(diss_tke./ro);
% ReL = (eke.*eke./eps).*ro./mu;
% ReLam = sqrt(abs(20*ReL/3));

I = 1:ni-1;
dxi = x(I+1,:)-x(I,:);
dxi(ni,:) = dxi(ni-1,:);
dyi = y(I+1,:)-y(I,:);
dyi(ni,:) = dyi(ni-1,:);
J = 1:nj-1;
dxj = x(:,J+1)-x(:,J);
dxj(:,nj) = dxj(:,nj-1);
dyj = y(:,J+1)-y(:,J);
dyj(:,nj) = dyj(:,nj-1);
dcell = sqrt(abs(dxi.*dyj - dxj.*dyi));
% J = dxi.*dyj - dxj.*dyi;


if(nb==1)
    [~,ii] = min( abs(mean(x,2) - xinlet) );
    iinlet = ii;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    Toav = (To(ii,jj+1)+To(ii,jj))*0.5;
    
    pav = (p(ii,jj+1)+p(ii,jj))*0.5;
    rusav = (rus(ii,jj+1)+rus(ii,jj))*0.5;
    ruvav = (ruv(ii,jj+1)+ruv(ii,jj))*0.5;
    roav = (ro(ii,jj+1)+ro(ii,jj))*0.5;
    
    ekeav = (eke(ii,jj+1)+eke(ii,jj))*0.5;

    poin = sum(poav.*dm);
    Toin = sum(Toav.*dm);
    roin = sum(roav.*dy);
    
    pin = sum(pav.*dy);
    rusin = sum(rusav.*dy);
    ruvin = sum(ruvav.*dy);
    ekein = sum(ekeav.*dy);
      
    mdot_in = sum(dm);
    YY = sum(dy);
    
    
    yinlet = y(ii,:)';
    ptinlet = po(ii,:)';
    ekeinlet = eke(ii,:)';
    
end
if(nb==2)
    %  [~,ii] = min( abs(mean(x,2) - xinlet) );
    ii = iinlet;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    Toav = (To(ii,jj+1)+To(ii,jj))*0.5;
    
    pav = (p(ii,jj+1)+p(ii,jj))*0.5;
    rusav = (rus(ii,jj+1)+rus(ii,jj))*0.5;
    ruvav = (ruv(ii,jj+1)+ruv(ii,jj))*0.5;
    roav = (ro(ii,jj+1)+ro(ii,jj))*0.5;

    ekeav = (eke(ii,jj+1)+eke(ii,jj))*0.5;
      
    poin = poin + sum(poav.*dm);
    Toin = Toin + sum(Toav.*dm);
    
    pin = pin + sum(pav.*dy);
    rusin = rusin + sum(rusav.*dy);
    ruvin = ruvin + sum(ruvav.*dy);
    roin = roin + sum(roav.*dy);
    ekein = ekein + sum(ekeav.*dy);
    
    mdot_in = mdot_in + sum(dm);
    YY = YY + sum(dy);
    
    poin = poin/mdot_in;
    Toin = Toin/mdot_in;
    
    pin = pin/YY;
    roin = roin/YY;
    vin = ruvin/mdot_in;
    uin = mdot_in/(YY*roin);
    ent_in = rusin/mdot_in;
    ekein = ekein/YY;
    
    yinlet = [yinlet; y(ii,:)'];
    ptinlet = [ptinlet; po(ii,:)'];
    ekeinlet = [ekeinlet; eke(ii,:)'];
    
end


if(nb==8)
[~,ii] = min( abs(mean(x,2) - xexit) );
    iexit = ii;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    roav = (ro(ii,jj+1)+ro(ii,jj))*0.5;    
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    pav = (p(ii,jj+1)+p(ii,jj))*0.5;
    
    rusav = (rus(ii,jj+1)+rus(ii,jj))*0.5;
    ruvav = (ruv(ii,jj+1)+ruv(ii,jj))*0.5;
    
    poex = sum(poav.*dm);
    pex = sum(pav.*dy);    
    rusex = sum(rusav.*dy);
    ruvex = sum(ruvav.*dy);
    roex =  sum(roav.*dy);
    mdot = sum(dm);
    dY = sum(dy);

end
if(nb==9)
%[~,ii] = min( abs(mean(x,2) - xexit) );
    ii = iexit;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    pav = (p(ii,jj+1)+p(ii,jj))*0.5;
        
    roav = (ro(ii,jj+1)+ro(ii,jj))*0.5;
    
    rusav = (rus(ii,jj+1)+rus(ii,jj))*0.5;
    ruvav = (ruv(ii,jj+1)+ruv(ii,jj))*0.5;

    poex = poex + sum(poav.*dm);
    pex = pex + sum(pav.*dy);
      
    roex = roex + sum(roav.*dy);
    
    rusex = rusex + sum(rusav.*dy);
    ruvex = ruvex + sum(ruvav.*dy);
    
    mdot = mdot + sum(dm);
    
    dY = dY + sum(dy);
    
    roex = roex/dY;
    poex = poex/mdot;
    pex = pex/dY;
    
    vex = ruvex/mdot;
    uex = mdot/(dY*roex);
    
end


if(nb==8)
[~,ii] = min( abs(mean(x,2) - xwake) );
    iwake = ii;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    mdot_wake = sum(dm);
    powake = sum(poav.*dm);     
end
if(nb==9)
%[~,ii] = min( abs(mean(x,2) - xwake) );
    ii = iwake;
    jj=1:nj-1;
    dy = y(ii,jj+1)-y(ii,jj);
    dm = (ru(ii,jj+1)+ru(ii,jj)).*dy*0.5;
    poav = (po(ii,jj+1)+po(ii,jj))*0.5;
    mdot_wake = mdot_wake + sum(dm);       
    powake = powake + sum(poav.*dm);
    powake = powake/mdot_wake;     
end



if(nb==8 || nb==9)
% figure(1)  
 ii = iwake;

if nb==8
ywake = y(ii,:)';
Twake = T(ii,:)';
ptwake = po(ii,:)';
uwake = u(ii,:)';
vwake = v(ii,:)';
swake = s(ii,:)';
pwake = p(ii,:)';
mwake = mach(ii,:)';
kwake = eke(ii,:)';
urmswake = urms(ii,:)';
vrmswake = vrms(ii,:)';
wrmswake = wrms(ii,:)';
else
ywake = [ywake; y(ii,:)'];
Twake = [Twake; T(ii,:)'];
ptwake = [ptwake; po(ii,:)'];
uwake = [uwake; u(ii,:)'];
vwake = [vwake; v(ii,:)'];
swake = [swake; s(ii,:)'];
pwake = [pwake; p(ii,:)'];
mwake = [mwake; mach(ii,:)'];
kwake = [kwake; eke(ii,:)'];
urmswake = [urmswake; urms(ii,:)'];
vrmswake = [vrmswake; vrms(ii,:)'];
wrmswake = [wrmswake; wrms(ii,:)'];
end


if(nb==8)
xwake8 = mean(x(ii,:));
elseif(nb==9)
xwake9 = mean(x(ii,:));
end

end



qbar.ro = ro;
qbar.ru = ru;
qbar.rv = rv;
qbar.rw = rw;
qbar.Et = Et;

qbar.ruu = ruu;
qbar.rvv = rvv;
qbar.rww = rww;

qbar.ruv = ruv;
qbar.ruw = ruw;
qbar.rvw = rvw; 

if(nb==3)

blparams.ISS = 0;
blparams.IPS = 0;
blparams.XX = 0;
blparams.YY = 0;
blparams.SS = 0;
blparams.x0 = 0;
blparams.y0 = 0;

blparams.d0 = 0;
blparams.d1 = 0;
blparams.d2 = 0;
blparams.d3 = 0;

blparams.d1i = 0;
blparams.d2i = 0;
blparams.d3i = 0;

blparams.d4 = 0;
blparams.d5 = 0;
%blparams.d6 = d6;
blparams.d7 = 0;
%blparams.d8 = d8;
%blparams.d9 = d9;
blparams.d10 = 0;
blparams.d11 = 0;
blparams.d12 = 0;
blparams.d13 = 0;
%blparams.d14 = d14;
blparams.d15 = 0;
blparams.d16 = 0;
blparams.d17 = 0;
blparams.ue = 0;
blparams.d3u3 = 0;
blparams.skinfriction = 0;
blparams.retheta = 0;    
blparams.tke_edge = 0; 
blparams.pwall = 0;
nsnow = 0;
end

if(nb == 3 || nb == 4 || nb == 5 || nb == 7)
blparams = integralparams(x,y,qbar,mu_ref,gam,rgas,vref,nb,nsnow,xLE,xTE,blparams);
end


flow{nb}.qbar = qbar;
% 

if(nb == 3 || nb == 4 || nb == 5 || nb == 7)
nsnow = nsnow + ni;
end
ifacemax = iface_tot + (nj-2)*(ni-1) + ni-1;
faces(iface_tot+1:ifacemax,1:4) = 0;

for j=1:nj-1
for i=1:ni-1
iface = iface_tot + (j-1)*(ni-1) + i;        
faces(iface,1) = nn + ni*(j-1) + i;
faces(iface,2) = nn + ni*j + i;
faces(iface,3) = nn + ni*j + i+1;
faces(iface,4) = nn + ni*(j-1) + i+1;       
end
end
% 
% 
for j=1:nj
for i=1:ni   
nn=nn+1;
if(nb == 3 || nb==4 || nb ==5 || nb == 7)
%ydist = y(i,j)-y(i,nj);
ydist = y(i,nj-2)-y(i,nj);
else
ydist = 0;
end
% % 
% % IPS(nn) = ydist < 0;
% % ISS(nn) = ydist > 0;
BLFLAG = 0;
PSFLAG = ydist < 0;
SSFLAG = ydist > 0;
if(PSFLAG || SSFLAG)
BLFLAG = j > blparams.edgej(i);
end
meanflow.X(nn) = x(i,j);
meanflow.Y(nn) = y(i,j);

meanflow.RO(nn) = ro(i,j);
meanflow.RU(nn) = ru(i,j);
meanflow.RV(nn) = rv(i,j);
meanflow.RW(nn) = rw(i,j);

meanflow.RUU(nn) = ruudash(i,j);
meanflow.RVV(nn) = rvvdash(i,j);
meanflow.RWW(nn) = rwwdash(i,j);
meanflow.RUV(nn) = ruvdash(i,j);
meanflow.P(nn) = p(i,j);

meanflow.TPROD(nn) = turb(i,j)/ro(i,j);
meanflow.EDGEWORK(nn) = turb2(i,j)/ro(i,j);
meanflow.TKEFLUX(nn) = div_tke(i,j);
meanflow.LOSS(nn) = (rvgrads(i,j) + rdsdt(i,j)).*T(i,j)./ro(i,j);
meanflow.LOSS_STEADY(nn) = diss_av(i,j)./ro(i,j);
meanflow.DISS(nn) = diss(i,j)/ro(i,j);
meanflow.KEFLUX(nn) = ke_flux(i,j);
meanflow.KEISENFLUX(nn) = ke_isen_flux(i,j);
meanflow.DIFFK(nn) = diffK(i,j)./ro(i,j);
meanflow.PS(nn) = PSFLAG*BLFLAG;
meanflow.SS(nn) = SSFLAG*BLFLAG;
meanflow.VORTZ(nn) = dudy(i,j) - dvdx(i,j);
meanflow.ENTFLUX(nn) = rvgrads(i,j);
meanflow.ENTPROD(nn) = diss(i,j)/T(i,j);

end
end
% 
iface_tot = iface_tot + (ni-1)*(nj-1);

end


save([dir,'faces'],'faces')

Tex = pex/(rgas*roex);
mu_ex = (mu_ref)*( ( 273.0 + 110.4 )./( Tex + 110.4 ) ).*((Tex/273).^1.5) ;
vel_ex =  sqrt(uex*uex + vex*vex);
Rey_exit = vel_ex*roex*cax/mu_ex;
M_exit = vel_ex/sqrt(gam*rgas*Tex);

Tex_isen = Toin*((pex/poin).^((gam-1)/gam));
vex_isen = sqrt((Toin - Tex_isen)*2*cp);

save([dir,'meanflow'],'meanflow','Tex_isen','vex_isen','mdot');

% suction surface dissipation
ISS = blparams.ISS == 1;
XS = blparams.XX(ISS);
YS = blparams.YY(ISS);
X0 = blparams.x0(ISS);
Y0 = blparams.y0(ISS);
D0 = blparams.d0(ISS);
UE = blparams.ue(ISS);
D1 = blparams.d1(ISS);
D2 = blparams.d2(ISS);
D3 = blparams.d3(ISS);

D1i = blparams.d1i(ISS);
D2i = blparams.d2i(ISS);
D3i = blparams.d3i(ISS);

RET = blparams.retheta(ISS);
UK = blparams.d4(ISS);
D5 = blparams.d5(ISS);
% D6 = blparams.d6(ISS);
% D9 = blparams.d9(ISS);
D10 = blparams.d10(ISS);
D11 = blparams.d11(ISS);
% D12 = blparams.d12(ISS);
D13 = blparams.d13(ISS);
% D14 = blparams.d14(ISS);
% D16 = blparams.d16(ISS);
% D17 = blparams.d17(ISS);
% D19 = blparams.d19(ISS);
% D20 = blparams.d20(ISS);
% D21 = blparams.d21(ISS);
% D22 = blparams.d22(ISS);
% D23 = blparams.d23(ISS);
% D24 = blparams.d24(ISS);
%VOL = volbl(ISS);
CF = blparams.skinfriction(ISS);
TUE = sqrt((blparams.tke_edge(ISS))*2/3)./UE;

% sort:
[XS,I] = sort(XS);
YS = YS(I);
X0 = X0(I);
Y0 = Y0(I);
D0 = D0(I);
UE = UE(I);

D1 = D1(I);
D2 = D2(I);
D3 = D3(I);

D1i = D1i(I);
D2i = D2i(I);
D3i = D3i(I);

RET = RET(I);
UK = UK(I);
D5 = D5(I);
% D6 = D6(I);
% D9 = D9(I);
D10 = D10(I);
D11 = D11(I);
% D12 = D12(I);
D13 = D13(I);
% D14 = D14(I);
% D16 = D16(I);
% D17 = D17(I);
% D19 = D19(I);
% D20 = D20(I);
% D21 = D21(I);
% D22 = D22(I);
% D23 = D23(I);
% D24 = D24(I);
%VOL = VOL(I);
CF = CF(I);
TUE = TUE(I);

SS(1) = 0;

N = length(XS);

% surface distance
for i = 2:N
dx = XS(i)-XS(i-1);
dy = YS(i)-YS(i-1);
ds = sqrt(dx*dx + dy*dy);
SS(i) = SS(i-1) + ds;
end

%dukds = dukds./((UE.^3))
%Cd =  D12./((UE.^3));
dukds = D11./((UE.^3));
D5 = D5./((UE.^3));
%D6 = D6./((UE.^3)); % entropy flux
%D9 = D9./((UE.^3)); % thermal term
D10 = D10./((UE.^3)); % Turbulent production term (apparent dissipation due to Rey stresses)
D13 = D13./((UE.^3)); % Dissipation computed from time-average strain
%D14 = D14./((UE.^3));
% 
% D16 = D16./((UE.^3));
% D17 = D17./((UE.^3));
% D19 = D19./((UE.^3)); % Turbulent work on boundaries- i.e. acceleration or deceleration of base flow by turbulence
% D20 = D20./((UE.^3)); % difference between production and turbulent work on boundaries
% D21 = D21./((UE.^3));
% D22 = D22./((UE.^3));
% D24 = D24./((UE.^3));


suction.xs = XS;
suction.ss = SS;
suction.d0= D0;
suction.d1 = D1;
suction.d2 = D2;
suction.d3 = D3;

suction.d1i = D1i;
suction.d2i = D2i;
suction.d3i = D3i;

suction.d5 = D5;
suction.d10 = D10;
suction.d13 = D13;
suction.dukds= dukds;
suction.ue = UE;
suction.tue = TUE;
suction.cf = CF;
suction.ret = RET;

exit.rex = Rey_exit;
exit.M = M_exit;
exit.vel = vel_ex;
exit.T = Tex;
exit.p = pex;
exit.visen = vex_isen;
exit.Tisen = Tex_isen;

exit.mdot = mdot;
exit.po = poex;
exit.u = uex;
exit.v = vex;
exit.s = rusex./mdot;

wake.y = ywake;
wake.T = Twake;
wake.po = ptwake;
wake.u = uwake;
wake.v = vwake;
wake.s = swake;
wake.p = pwake;
wake.m = mwake;
wake.urms = urmswake;
wake.vrms = vrmswake;
wake.wrms = wrmswake;
wake.tke = kwake;

inlet.mdot = mdot_in;
inlet.po = poin;
inlet.To = Toin;
inlet.p = pin;
inlet.ro = roin;
inlet.u = uin;
inlet.v = vin;
inlet.tke = ekein;
inlet.s = ent_in;

inlet.y = yinlet;
inlet.Pt = ptinlet;
inlet.Tke = ekeinlet;

blparams.cp = (poin - blparams.pwall)/(poin - pex);

return


