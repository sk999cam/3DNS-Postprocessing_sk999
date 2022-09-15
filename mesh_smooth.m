function [blk] = mesh_smooth(blk,next_block,next_patch,corner,pitch,msmooths,x_prof,y_prof,ywall,iupdate)

if nargin < 10
    iupdate = 10;
end


NB = length(blk);
for ib=1:NB
    [nib{ib},njb{ib}] = size(blk{ib}.x);
    NI{ib} = nib{ib};
    NJ{ib} = njb{ib};
end

% ensure wall profile remains fixed;
s_prof = curve_length(x_prof,y_prof);
s_len = s_prof(end);

s_prof2 = [s_prof-s_len s_prof(2:end-1) s_prof+s_len];
x_prof2 = [x_prof x_prof(2:end-1) x_prof];
y_prof2 = [y_prof y_prof(2:end-1) y_prof];

for mm=1:msmooths
    %
    % % fix periodics
    % for ib=1:NB
    % [blk]=make_periodic_smooth(blk,next_block,next_patch,ib,pitch);
    % %[blk]=make_periodic(blk,next_block,next_patch,ib,pitch);
    % end

    if mod(mm,iupdate) == 0
        fprintf('Iteration %d\n',mm)
    end

    % call multiblock to find the block boundary interfaces
    ng = 2;
    for ib=1:NB
        [up{ib},dn{ib}] = multiblock(blk,next_block,next_patch,ib,pitch,ng);
    end

    % smooth routine


    for ib=1:NB
        im_next_block = next_block{ib}.im;
        ip_next_block = next_block{ib}.ip;
        jm_next_block = next_block{ib}.jm;
        jp_next_block = next_block{ib}.jp;
        
        im_next_patch = next_patch{ib}.im;
        ip_next_patch = next_patch{ib}.ip;
        jm_next_patch = next_patch{ib}.jm;
        jp_next_patch = next_patch{ib}.jp;
        
        xnew = blk{ib}.x ;
        ynew = blk{ib}.y ;
        
        im_wall = im_next_block ==  0 && im_next_patch == 3;
        ip_wall = ip_next_block ==  0 && ip_next_patch == 3;
        jm_wall = jm_next_block ==  0 && jm_next_patch == 3;
        jp_wall = jp_next_block ==  0 && jp_next_patch == 3;
        
        
        im_in = im_next_block ==  0 && im_next_patch == 1;
        ip_ex = ip_next_block ==  0 && ip_next_patch == 2;
        
        
        [ni_new,nj_new] = size(xnew);
        
        xdum = xnew;
        ydum = ynew;
        x = xnew;
        y = ynew;
        
        is = 1;
        ie = ni_new;
        js = 1;
        je = nj_new;
        
        
        if(im_next_block==0)% || im_wall)
            is = 2;
        end
        
        if(jm_next_block==0)% || jm_wall)
            js = 2;
        end
        
        if(ip_next_block==0)%  || ip_wall)
            ie = ni_new-1;
        end
        
        if(jp_next_block==0)% || jp_wall)
            je = nj_new-1;%nj_new-1;
        end
        
        
        
        for j=js:je
            for i=is:ie
            
                if(i==1)
                    xim = up{ib}.xi(ng,j);
                    yim = up{ib}.yi(ng,j);
                else
                    xim = x(i-1,j);
                    yim = y(i-1,j);
                end
            
                if(i==ni_new)
                    xip = dn{ib}.xi(ng,j);
                    yip = dn{ib}.yi(ng,j);
                else
                    xip = x(i+1,j);
                    yip = y(i+1,j);
                end
            
            
                if(j==1)
                    xjm = up{ib}.xj(i,ng);
                    yjm = up{ib}.yj(i,ng);
                else
                    xjm = x(i,j-1);
                    yjm = y(i,j-1);
                end
            
                if(j==nj_new)
                    xjp = dn{ib}.xj(i,ng);
                    yjp = dn{ib}.yj(i,ng);
                else
                    xjp = x(i,j+1);
                    yjp = y(i,j+1);
                end
            
                dxi(i,j) = (-xim + xip);
                dyi(i,j) = (-yim + yip);
            
                dxj(i,j) = (-xjm + xjp);
                dyj(i,j) = (-yjm + yjp);
            
                %     if(jp_wall && j==nj_new-1)
                %     dxj(i,j) = dxj(i,j)*ywall/sqrt(dxj(i,j)*dxj(i,j) + dyj(i,j)*dyj(i,j));
                %     dyj(i,j) = dyj(i,j)*ywall/sqrt(dxj(i,j)*dxj(i,j) + dyj(i,j)*dyj(i,j));
                %     end
            
            
                gamma(i,j) = dxi(i,j)*dxi(i,j) + dyi(i,j)*dyi(i,j);
                alpha(i,j) = dxj(i,j)*dxj(i,j) + dyj(i,j)*dyj(i,j);
                beta(i,j)  = dxj(i,j)*dxi(i,j) + dyj(i,j)*dyi(i,j);
            
            
            
                if(i==1 && j>1 && j<nj_new)
                    d2xij =  (x(i+1,j+1) - x(i+1,j-1) + up{ib}.xi(ng,j-1) - up{ib}.xi(ng,j+1));
                    d2yij =  (y(i+1,j+1) - y(i+1,j-1) + up{ib}.yi(ng,j-1) - up{ib}.yi(ng,j+1));
                elseif(i==ni_new && j>1 && j<nj_new)
                    d2xij =  (dn{ib}.xi(ng,j+1) - dn{ib}.xi(ng,j-1) + x(i-1,j-1) - x(i-1,j+1));
                    d2yij =  (dn{ib}.yi(ng,j+1) - dn{ib}.yi(ng,j-1) + y(i-1,j-1) - y(i-1,j+1));
                elseif(j==1 && i>1 && i<ni_new)
                    d2xij =  (x(i+1,j+1) - up{ib}.xj(i+1,ng) + up{ib}.xj(i-1,ng) - x(i-1,j+1));
                    d2yij =  (y(i+1,j+1) - up{ib}.yj(i+1,ng) + up{ib}.yj(i-1,ng) - y(i-1,j+1));
                elseif(j==nj_new && i>1 && i<ni_new)
                    d2xij =  (dn{ib}.xj(i+1,ng) - x(i+1,j-1) + x(i-1,j-1) - dn{ib}.xj(i-1,ng));
                    d2yij =  (dn{ib}.yj(i+1,ng) - y(i+1,j-1) + y(i-1,j-1) - dn{ib}.yj(i-1,ng));
                elseif( (i>1 && i<ni_new) && (j>1 && j<nj_new))
                    d2xij =  (x(i+1,j+1) - x(i+1,j-1) + x(i-1,j-1) - x(i-1,j+1));
                    d2yij =  (y(i+1,j+1) - y(i+1,j-1) + y(i-1,j-1) - y(i-1,j+1));
                else
                    d2xij =  0;
                    d2yij =  0;
                end
            
            
            
                xnew(i,j) = 0.5*((xim + xip)*alpha(i,j) - beta(i,j)*d2xij*0.5 ...
                    +(xjm + xjp)*gamma(i,j) )...
                    /(alpha(i,j) + gamma(i,j)) ;
                ynew(i,j) = 0.5*((yim + yip)*alpha(i,j) - beta(i,j)*d2yij*0.5 ...
                    +(yjm + yjp)*gamma(i,j) )...
                    /(alpha(i,j) + gamma(i,j)) ;
            
            
            end
        end


        blk{ib}.x = xnew;
        blk{ib}.y = ynew;


    
    
        
        
        
        % Apply boundary conditions
        
        % build o-grid
        if(jp_wall)
            
            xprof = x(:,end);
            yprof = y(:,end);
            sprof = curve_length(xprof,yprof);
            s_len_local = sprof(end) - sprof(1);

            si = linspace(sprof(1)-s_len_local*0.1,sprof(ni_new)+s_len_local*0.1,ni_new*100);
            [sprof,iu] = unique(sprof);
            xprof = xprof(iu);
            yprof = yprof(iu);
            
            % find start point on original profile
            [~,istart]=min( sqrt( (x_prof-xprof(1)).^2 + (y_prof-yprof(1)).^2 ) );
            if ((ib==4 || ib==7) && (NB == 9))
                 [~,istart]=min( sqrt( (x_prof-xprof(end)).^2 + (y_prof-yprof(end)).^2 ) );
            end
            
            if ((ib==5 || ib==9) && (NB == 12))
                 [~,istart]=min( sqrt( (x_prof-xprof(end)).^2 + (y_prof-yprof(end)).^2 ) );
            end

            
            sstart = s_prof(istart);% + sqrt((x_prof(istart)-xprof(1))^2 + (y_prof(istart)-yprof(1))^2);
            si = si + sstart;
            
            si(si>s_len) = si(si>s_len)-s_len;
            si(si<0) = si(si<0)+s_len;
            
            xi = interp1(s_prof2,x_prof2,si,'spline');
            yi = interp1(s_prof2,y_prof2,si,'spline');
            
            for i=1:ni_new
                % find nearest wall point to set orthogonality
                d = sqrt((x(i,nj_new-1)-xi).^2 + (y(i,nj_new-1)-yi).^2);
                [~,inorm]=min(d);
                xwall_now(i) = xi(inorm);
                ywall_now(i) = yi(inorm);
                xnew(i,nj_new) = xi(inorm);
                ynew(i,nj_new) = yi(inorm);
            end
            
            %
            % now drive near wall distance to ywall
            for i=1:ni_new
                ynorm = curve_length(xnew(i,:),ynew(i,:));
                % get expansion factor and spacing
                fex = fexpan(ynorm(nj_new)/ywall,nj_new);
                fy = spacing(nj_new,1/fex,0);
                
                % fit spline to find new wall normal line
                ynorm = ynorm/ynorm(nj_new);
                yni = fy;
                xnew(i,:) = interp1(ynorm,xnew(i,:),yni,'spline');
                ynew(i,:) = interp1(ynorm,ynew(i,:),yni,'spline');
            end

        end

        
        
        if(im_in)
            ynew(1,:) = ynew(2,:);
        end
        
        if(ip_ex)
            ynew(end,:) = ynew(end-1,:);
        end
        
        
        
        
        blk{ib}.x = xnew;
        blk{ib}.y = ynew;

        
        
        clear x y

    end


    %
    % make periodic
    for ib=1:NB
        [blk]=make_periodic(blk,next_block,next_patch,ib,pitch);
    end

    %
    % fix periodics
    for ib=1:NB
        [blk]=make_periodic_smooth2(blk,next_block,next_patch,ib,pitch);
    end

    %
    % Now corner treatment
    % in this section the grid points where i=j are shrunk toward the corner
    % point
    ncorner = length(corner);
    for n=1:ncorner

        xcor{n} = 0.0;
        ycor{n} = 0.0;
        dely = 0.0;

        for m=1:corner{n}.Nb
            ib = corner{n}.block{m};
            ic = corner{n}.i{m};
            jc = corner{n}.j{m};

            if ic~=1; ic = NI{ib}; end
            if jc~=1; jc = NJ{ib}; end

            if(m==1)
                yc0 = blk{ib}.y(ic,jc);
            end

            yoffset = 0;
            if( (blk{ib}.y(ic,jc) - yc0) > pitch*0.1); yoffset = -pitch; end
            if( (blk{ib}.y(ic,jc) - yc0) <-pitch*0.1); yoffset =  pitch; end

            xcor{n} = xcor{n} +  blk{ib}.x(ic,jc);
            ycor{n} = ycor{n} +  blk{ib}.y(ic,jc) + yoffset;

        end

        xc{n} = xcor{n}/corner{n}.Nb;
        yc{n} = ycor{n}/corner{n}.Nb;
        yc_0{n} = yc0;

        xx = xc{n};
        yy = yc{n};

        nic =1e12;
        for m=1:corner{n}.Nb
            ib = corner{n}.block{m};
            nic_now = min([NI{ib}-1,NJ{ib}-1]);
            if(nic_now<nic)
                nic = nic_now;
            end

        end

        nic = floor(nic/2);

        for m=1:corner{n}.Nb
            ib = corner{n}.block{m};
            ic = corner{n}.i{m};
            jc = corner{n}.j{m};

            if ic~=1; ic = NI{ib}; end
            if jc~=1; jc = NJ{ib}; end

            yoffset = 0;
            if( (blk{ib}.y(ic,jc) - yc_0{n}) > pitch*0.1); yoffset = -pitch; end
            if( (blk{ib}.y(ic,jc) - yc_0{n}) <-pitch*0.1); yoffset =  pitch; end

            xx = xc{n};
            yy = yc{n}-yoffset;

            [ni_new,nj_new]=size(blk{ib}.x);

            blk{ib}.x(ic,jc) = xx;
            blk{ib}.y(ic,jc) = yy;
            %
            nic=4;
            fc=ones(nic,nic);
            if(corner{n}.Nb==5)
                fc = fc*nic;
                % for id=1:nic
                %     fc(id,id) = id;
                % end
                fc(1:nic,1)=1:nic;
                fc(1,1:nic)=1:nic;
                fc=fc/nic;
                fc = fc.^0.25;%fc.^0.25; % change the power to change shrinkage to corner point
            end

            ii=1:nic;
            jj=1:nic;
            if(ic~=1); ii=ni_new:-1:ni_new-nic+1; end
            if(jc~=1); jj=nj_new:-1:nj_new-nic+1; end

            blk{ib}.x(ii,jj) = xx + (blk{ib}.x(ii,jj)-xx).*fc;
            blk{ib}.y(ii,jj) = yy + (blk{ib}.y(ii,jj)-yy).*fc;


        end

    end

    %
    % make periodic
    for ib=1:NB
        [blk]=make_periodic(blk,next_block,next_patch,ib,pitch);
    end

end



return