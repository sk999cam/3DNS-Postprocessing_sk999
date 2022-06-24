function [blk] = make_periodic_smooth2(blk,next_block,next_patch,nbnow,pitch)

ng = 1;

up.xi = 0;
up.yi = 0;
dn.xi = 0;
dn.yi = 0;

up.xj = 0;
up.yj = 0;
dn.xj = 0;
dn.yj = 0;



ib = nbnow;

im_next_block = next_block{ib}.im;
ip_next_block = next_block{ib}.ip;
jm_next_block = next_block{ib}.jm;
jp_next_block = next_block{ib}.jp;

im_next_patch = next_patch{ib}.im;
ip_next_patch = next_patch{ib}.ip;
jm_next_patch = next_patch{ib}.jm;
jp_next_patch = next_patch{ib}.jp;

f1 = 0.5; %16.0/30.0
f2 = 0.0; %-1.0/30.0
f0 = 0.95;
f1 = f1*(1.0-f0);
f2 = f2*(1.0-f0);

if(im_next_block~=0)

    if(im_next_patch==1)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end

        xup = blk{im_next_block}.x(2,:);
        yup = blk{im_next_block}.y(2,:)+yoffset;
        xup2 = blk{im_next_block}.x(3,:);
        yup2 = blk{im_next_block}.y(3,:)+yoffset;

        xnow = blk{ib}.x(1,:);
        ynow = blk{ib}.y(1,:);

        blk{ib}.x(1,:) = xnow*f0 + ...
            (xup + blk{ib}.x(2,:))*f1+(xup2 + blk{ib}.x(3,:))*f2;
        blk{ib}.y(1,:) = ynow*f0 + ...
            (yup + blk{ib}.y(2,:))*f1+(yup2 + blk{ib}.y(3,:))*f2;

        blk{im_next_block}.x(1,:) = blk{ib}.x(1,:);
        blk{im_next_block}.y(1,:) = blk{ib}.y(1,:)-yoffset;
    end

    if(im_next_patch==2)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end

        xup = blk{im_next_block}.x(end-1,:);
        yup = blk{im_next_block}.y(end-1,:)+yoffset;
        xup2 = blk{im_next_block}.x(end-2,:);
        yup2 = blk{im_next_block}.y(end-2,:)+yoffset;

        xnow = blk{ib}.x(1,:);
        ynow = blk{ib}.y(1,:);

        blk{ib}.x(1,:) = xnow*f0 + ...
            (xup + blk{ib}.x(2,:))*f1+(xup2 + blk{ib}.x(3,:))*f2;
        blk{ib}.y(1,:) = ynow*f0 + ...
            (yup + blk{ib}.y(2,:))*f1+(yup2 + blk{ib}.y(3,:))*f2;

        blk{im_next_block}.x(end,:) = blk{ib}.x(1,:);
        blk{im_next_block}.y(end,:) = blk{ib}.y(1,:)-yoffset;
    end

    if(im_next_patch==3)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)') >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)') < -0.1*pitch ); yoffset =-pitch; end

        xup = blk{im_next_block}.x(:,2)';
        yup = blk{im_next_block}.y(:,2)'+yoffset;
        xup2 = blk{im_next_block}.x(:,3)';
        yup2 = blk{im_next_block}.y(:,3)'+yoffset;

        xnow = blk{ib}.x(1,:);
        ynow = blk{ib}.y(1,:);

        blk{ib}.x(1,:) = xnow*f0 + ...
            (xup + blk{ib}.x(2,:))*f1+(xup2 + blk{ib}.x(3,:))*f2;
        blk{ib}.y(1,:) = ynow*f0 + ...
            (yup + blk{ib}.y(2,:))*f1+(yup2 + blk{ib}.y(3,:))*f2;

        blk{im_next_block}.x(:,1) = blk{ib}.x(1,:)';
        blk{im_next_block}.y(:,1) = blk{ib}.y(1,:)'-yoffset;
    end

    if(im_next_patch==4)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,end)') >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,end)') < -0.1*pitch ); yoffset =-pitch; end

        xup = blk{im_next_block}.x(:,end-1)';
        yup = blk{im_next_block}.y(:,end-1)'+yoffset;
        xup2 = blk{im_next_block}.x(:,end-2)';
        yup2 = blk{im_next_block}.y(:,end-2)'+yoffset;

        xnow = blk{ib}.x(1,:);
        ynow = blk{ib}.y(1,:);

        blk{ib}.x(1,:) = xnow*f0 + ...
            (xup + blk{ib}.x(2,:))*f1+(xup2 + blk{ib}.x(3,:))*f2;
        blk{ib}.y(1,:) = ynow*f0 + ...
            (yup + blk{ib}.y(2,:))*f1+(yup2 + blk{ib}.y(3,:))*f2;

        blk{im_next_block}.x(:,end) = blk{ib}.x(1,:)';
        blk{im_next_block}.y(:,end) = blk{ib}.y(1,:)'-yoffset;
    end

    up.xi = xup;
    up.yi = yup;

end

if(ip_next_block~=0)

    if(ip_next_patch==1)
        yoffset = 0.0;
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{ip_next_block}.x(2,:);
        yup = blk{ip_next_block}.y(2,:)+yoffset;
        xup2 = blk{ip_next_block}.x(3,:);
        yup2 = blk{ip_next_block}.y(3,:)+yoffset;

        xnow = blk{ib}.x(end,:);
        ynow = blk{ib}.y(end,:);

        blk{ib}.x(end,:) = xnow*f0 + ...
            (xup + blk{ib}.x(end-1,:))*f1+(xup2 + blk{ib}.x(end-2,:))*f2;
        blk{ib}.y(end,:) = ynow*f0 + ...
            (yup + blk{ib}.y(end-1,:))*f1+(yup2 + blk{ib}.y(end-2,:))*f2;

        blk{ip_next_block}.x(1,:) = blk{ib}.x(end,:);
        blk{ip_next_block}.y(1,:) = blk{ib}.y(end,:)-yoffset;
    end

    if(ip_next_patch==2)
        yoffset = 0.0;
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{ip_next_block}.x(end-1,:);
        yup = blk{ip_next_block}.y(end-1,:)+yoffset;
        xup2 = blk{ip_next_block}.x(end-2,:);
        yup2 = blk{ip_next_block}.y(end-2,:)+yoffset;

        xnow = blk{ib}.x(end,:);
        ynow = blk{ib}.y(end,:);

        blk{ib}.x(end,:) = xnow*f0 + ...
            (xup + blk{ib}.x(end-1,:))*f1+(xup2 + blk{ib}.x(end-2,:))*f2;
        blk{ib}.y(end,:) = ynow*f0 + ...
            (yup + blk{ib}.y(end-1,:))*f1+(yup2 + blk{ib}.y(end-2,:))*f2;

        blk{ip_next_block}.x(end,:) = blk{ib}.x(end,:);
        blk{ip_next_block}.y(end,:) = blk{ib}.y(end,:)-yoffset;
    end

    if(ip_next_patch==3)
        yoffset = 0.0;
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)') >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)') < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{ip_next_block}.x(:,2)';
        yup = blk{ip_next_block}.y(:,2)'+yoffset;
        xup2 = blk{ip_next_block}.x(:,3)';
        yup2 = blk{ip_next_block}.y(:,3)'+yoffset;

        xnow = blk{ib}.x(end,:);
        ynow = blk{ib}.y(end,:);

        blk{ib}.x(end,:) = xnow*f0 + ...
            (xup + blk{ib}.x(end-1,:))*f1+(xup2 + blk{ib}.x(end-2,:))*f2;
        blk{ib}.y(end,:) = ynow*f0 + ...
            (yup + blk{ib}.y(end-1,:))*f1+(yup2 + blk{ib}.y(end-2,:))*f2;

        blk{ip_next_block}.x(:,1) = blk{ib}.x(end,:)';
        blk{ip_next_block}.y(:,1) = blk{ib}.y(end,:)'-yoffset;
    end

    if(ip_next_patch==4)
        yoffset = 0.0;
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,end)') >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,end)') < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{ip_next_block}.x(:,end-1)';
        yup = blk{ip_next_block}.y(:,end-1)'+yoffset;
        xup2 = blk{ip_next_block}.x(:,end-2)';
        yup2 = blk{ip_next_block}.y(:,end-2)'+yoffset;

        xnow = blk{ib}.x(end,:);
        ynow = blk{ib}.y(end,:);

        blk{ib}.x(end,:) = xnow*f0 + ...
            (xup + blk{ib}.x(end-1,:))*f1+(xup2 + blk{ib}.x(end-2,:))*f2;
        blk{ib}.y(end,:) = ynow*f0 + ...
            (yup + blk{ib}.y(end-1,:))*f1+(yup2 + blk{ib}.y(end-2,:))*f2;

        blk{ip_next_block}.x(:,end) = blk{ib}.x(end,:)';
        blk{ip_next_block}.y(:,end) = blk{ib}.y(end,:)'-yoffset;
    end

    dn.xi = xup;
    dn.yi = yup;

end




if(jm_next_block~=0)

    if(jm_next_patch==1)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jm_next_block}.x(2,:)';
        yup = blk{jm_next_block}.y(2,:)'+yoffset;
        xup2 = blk{jm_next_block}.x(3,:)';
        yup2 = blk{jm_next_block}.y(3,:)'+yoffset;

        xnow = blk{ib}.x(:,1);
        ynow = blk{ib}.y(:,1);

        blk{ib}.x(:,1) = xnow*f0 + ...
            (xup + blk{ib}.x(:,2))*f1+(xup2 + blk{ib}.x(:,3))*f2;
        blk{ib}.y(:,1) = ynow*f0 + ...
            (yup + blk{ib}.y(:,2))*f1+(yup2 + blk{ib}.y(:,3))*f2;

        blk{jm_next_block}.x(1,:) = blk{ib}.x(:,1)';
        blk{jm_next_block}.y(1,:) = blk{ib}.y(:,1)'-yoffset;
    end

    if(jm_next_patch==2)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jm_next_block}.x(end-1,:)';
        yup = blk{jm_next_block}.y(end-1,:)'+yoffset;
        xup2 = blk{jm_next_block}.x(end-2,:)';
        yup2 = blk{jm_next_block}.y(end-2,:)'+yoffset;

        xnow = blk{ib}.x(:,1);
        ynow = blk{ib}.y(:,1);

        blk{ib}.x(:,1) = xnow*f0 + ...
            (xup + blk{ib}.x(:,2))*f1+(xup2 + blk{ib}.x(:,3))*f2;
        blk{ib}.y(:,1) = ynow*f0 + ...
            (yup + blk{ib}.y(:,2))*f1+(yup2 + blk{ib}.y(:,3))*f2;

        blk{jm_next_block}.x(end,:) = blk{ib}.x(:,1)';
        blk{jm_next_block}.y(end,:) = blk{ib}.y(:,1)'-yoffset;
    end

    if(jm_next_patch==3)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jm_next_block}.x(:,2);
        yup = blk{jm_next_block}.y(:,2)+yoffset;
        xup2 = blk{jm_next_block}.x(:,3);
        yup2 = blk{jm_next_block}.y(:,3)+yoffset;

        xnow = blk{ib}.x(:,1);
        ynow = blk{ib}.y(:,1);

        blk{ib}.x(:,1) = xnow*f0 + ...
            (xup + blk{ib}.x(:,2))*f1+(xup2 + blk{ib}.x(:,3))*f2;
        blk{ib}.y(:,1) = ynow*f0 + ...
            (yup + blk{ib}.y(:,2))*f1+(yup2 + blk{ib}.y(:,3))*f2;

        blk{jm_next_block}.x(:,1) = blk{ib}.x(:,1);
        blk{jm_next_block}.y(:,1) = blk{ib}.y(:,1)-yoffset;
    end

    if(jm_next_patch==4)
        yoffset = 0.0;
        if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,end)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,end)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jm_next_block}.x(:,end-1);
        yup = blk{jm_next_block}.y(:,end-1)+yoffset;
        xup2 = blk{jm_next_block}.x(:,end-2);
        yup2 = blk{jm_next_block}.y(:,end-2)+yoffset;

        xnow = blk{ib}.x(:,1);
        ynow = blk{ib}.y(:,1);

        blk{ib}.x(:,1) = xnow*f0 + ...
            (xup + blk{ib}.x(:,2))*f1+(xup2 + blk{ib}.x(:,3))*f2;
        blk{ib}.y(:,1) = ynow*f0 + ...
            (yup + blk{ib}.y(:,2))*f1+(yup2 + blk{ib}.y(:,3))*f2;

        blk{jm_next_block}.x(:,end) = blk{ib}.x(:,1);
        blk{jm_next_block}.y(:,end) = blk{ib}.y(:,1)-yoffset;
    end

    up.xj = xup;
    up.yj = yup;

end


if(jp_next_block~=0)

    if(jp_next_patch==1)
        yoffset = 0.0;
        if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jp_next_block}.x(2,:)';
        yup = blk{jp_next_block}.y(2,:)'+yoffset;
        xup2 = blk{jp_next_block}.x(3,:)';
        yup2 = blk{jp_next_block}.y(3,:)'+yoffset;

        xnow = blk{ib}.x(:,end);
        ynow = blk{ib}.y(:,end);

        blk{ib}.x(:,end) = xnow*f0 + ...
            (xup + blk{ib}.x(:,end-1))*f1+(xup2 + blk{ib}.x(:,end-2))*f2;
        blk{ib}.y(:,end) = ynow*f0 + ...
            (yup + blk{ib}.y(:,end-1))*f1+(yup2 + blk{ib}.y(:,end-2))*f2;

        blk{jp_next_block}.x(1,:) = blk{ib}.x(:,end)';
        blk{jp_next_block}.y(1,:) = blk{ib}.y(:,end)'-yoffset;
    end

    if(jp_next_patch==2)
        yoffset = 0.0;
        if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jp_next_block}.x(end-1,:)';
        yup = blk{jp_next_block}.y(end-1,:)'+yoffset;
        xup2 = blk{jp_next_block}.x(end-2,:)';
        yup2 = blk{jp_next_block}.y(end-2,:)'+yoffset;

        xnow = blk{ib}.x(:,end);
        ynow = blk{ib}.y(:,end);

        blk{ib}.x(:,end) = xnow*f0 + ...
            (xup + blk{ib}.x(:,end-1))*f1+(xup2 + blk{ib}.x(:,end-2))*f2;
        blk{ib}.y(:,end) = ynow*f0 + ...
            (yup + blk{ib}.y(:,end-1))*f1+(yup2 + blk{ib}.y(:,end-2))*f2;

        blk{jp_next_block}.x(end,:) = blk{ib}.x(:,end)';
        blk{jp_next_block}.y(end,:) = blk{ib}.y(:,end)'-yoffset;
    end

    if(jp_next_patch==3)
        yoffset = 0.0;
        if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jp_next_block}.x(:,2);
        yup = blk{jp_next_block}.y(:,2)+yoffset;
        xup2 = blk{jp_next_block}.x(:,3);
        yup2 = blk{jp_next_block}.y(:,3)+yoffset;

        xnow = blk{ib}.x(:,end);
        ynow = blk{ib}.y(:,end);

        blk{ib}.x(:,end) = xnow*f0 + ...
            (xup + blk{ib}.x(:,end-1))*f1+(xup2 + blk{ib}.x(:,end-2))*f2;
        blk{ib}.y(:,end) = ynow*f0 + ...
            (yup + blk{ib}.y(:,end-1))*f1+(yup2 + blk{ib}.y(:,end-2))*f2;

        blk{jp_next_block}.x(:,1) = blk{ib}.x(:,end);
        blk{jp_next_block}.y(:,1) = blk{ib}.y(:,end)-yoffset;
    end

    if(jp_next_patch==4)
        yoffset = 0.0;
        if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,end)) >  0.1*pitch ); yoffset = pitch; end
        if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,end)) < -0.1*pitch ); yoffset =-pitch; end
        xup = blk{jp_next_block}.x(:,end-1);
        yup = blk{jp_next_block}.y(:,end-1)+yoffset;
        xup2 = blk{jp_next_block}.x(:,end-2);
        yup2 = blk{jp_next_block}.y(:,end-2)+yoffset;

        xnow = blk{ib}.x(:,end);
        ynow = blk{ib}.y(:,end);

        blk{ib}.x(:,end) = xnow*f0 + ...
            (xup + blk{ib}.x(:,end-1))*f1+(xup2 + blk{ib}.x(:,end-2))*f2;
        blk{ib}.y(:,end) = ynow*f0 + ...
            (yup + blk{ib}.y(:,end-1))*f1+(yup2 + blk{ib}.y(:,end-2))*f2;

        blk{jp_next_block}.x(:,end) = blk{ib}.x(:,end);
        blk{jp_next_block}.y(:,end) = blk{ib}.y(:,end)-yoffset;
    end

    dn.xj = xup;
    dn.yj = yup;

end






return