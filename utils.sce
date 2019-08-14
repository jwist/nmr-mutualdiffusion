function [Spectrum,GH] = R_data(path,GH),
    for i = 1:size(GH,2)
        idr=mopen(path+string(GH(1,i))+'/acqus','r');
        acqus = mgetl(idr,-1); mclose(idr);
        [w,r]=grep(acqus,'##$SPOFFS');
        [a,b,c,d]=regexp(acqus(w+1),'/(?P<digit>\D*\d+)\s+/');
        GH(2,i) = strtod(d(1));
        idr=mopen(path+string(GH(1,i))+'/pdata/1/1r','rb');
        [Spectrum(:,i),count_r] = mtlb_fread(idr, TD/2,'int32'); mclose(idr);
        idr = mopen(path+string(GH(1,i))+'/pdata/1/procs','r');
        procs = mgetl(idr,-1); mclose(idr);
        [w,r]=grep(procs,'##$NC_proc');
        [a,b,c,d]=regexp(procs(w),'/(?P<name>\w+)=\s+(?P<digit>\D*\d+)/');
        fa = strtod(d(2));
        Spectrum(:,i)= Spectrum(:,i).*(2^fa);
    end
endfunction

function GH=Order(GH,up_or_down),
    if up_or_down == 1 then
        for i = 1:size(GH,2)
            for j = 1:size(GH,2)-i
                if GH(2,j) > GH(2,j+1)
                    temp = GH(:,j);
                    GH(:,j) = GH(:,j+1);
                    GH(:,j+1) = temp;
                end
            end
        end
    elseif up_or_down == -1 then
        for i = 1:size(GH,2)
            for j = 1:size(GH,2)-i
                if GH(2,j) < GH(2,j+1)
                    temp = GH(:,j);
                    GH(:,j) = GH(:,j+1);
                    GH(:,j+1) = temp;
                end
            end
        end
    else
        error('the varible -up_or_down-  have not an avaliable value')
    end
    GH(3,:) = GH(1,:)-(folder_ini-1);
endfunction

function G_spectra(Spectrum,GH,delta_ppm,delta_int,style,up_or_down),
    a = gca();a.axes_reverse = ["on","off","off"] //scale as reference
    xtitle('','ppm','intensity');
    if (up_or_down == 1) | (up_or_down==-1) then
        GH = Order(GH,up_or_down); jump = 1;
        for i = GH(3,:)
            plot(ppm(1:$)'+delta_ppm*jump,Spectrum(1:$,i)+delta_int*jump,style);
            jump = jump + 1;
        end
    elseif up_or_down == 0 then
        for i = 1:size(Spectrum,2)
            plot(ppm(1:$)'+delta_ppm*(i-1),Spectrum(1:$,i)+delta_int*(i-1),style);
        end
    else
        error('the boolean varible -rand_spoffs-  have not an avaliable value')
    end
endfunction

function integral = R_integral(path,GH)
    counter = 1;
    for i = 1:size(GH,2)
        idr=mopen(path+string(GH(1,i))+'/pdata/1/intrng','r');
        intrng_temp = mgetl(idr,-1);mclose(idr);
        for j = 2:size(intrng_temp,1)
            [a,b,c,d]=regexp(intrng_temp(j),'/(?P<low_field_limit>\D*\d+\D*\d+)\s+(?P<high_field_limit>\D*\d+\D*\d+)/');
            intrng(1)=(find(ppm<strtod(d(1)),1))-1;
            intrng(2)=(find(ppm<strtod(d(2)),1))-1;
            I=sum(Spectrum(intrng(1):intrng(2),i));
            integral(counter,:) = [I,strtod(d(1)),strtod(d(2)),intrng(1),intrng(2),GH(2,i),GH(1,i)];
            counter = counter +1;
        end
    end
endfunction

function G_integral(list_intrng,integral)
    for l =1:size(list_intrng,1)
        for i = 1:size(integral,1)
            if integral(i,2) < list_intrng(l,1) then
                if integral(i,3) > list_intrng(l,2) then
                    plot2d(integral(i,6),integral(i,1),style=-list_intrng(l,3))
                end
            end
        end
    end
endfunction

function integral = Export_integral(list_intrng,integral)
    for j =1:size(list_intrng,1)
        for i = 1:size(integral,1)
            if integral(i,2) < list_intrng(j,1) then
                if integral(i,3) > list_intrng(j,2) then
                    integral(j,i,:) = [integral(i,1),integral(i,6),integral(i,7)];
                end
            end
        end
    end
endfunction

function integral = Integrate(list_intrng,Spectrum,GH,ppm)
    GH = Order(GH,1)
    for j = 1:size(list_intrng,1)
        for i = 1:size(Spectrum,2)
            intrng(1)=(find(ppm<list_intrng(j,1),1))-1;
            intrng(2)=(find(ppm<list_intrng(j,2),1))-1;
            I=sum(Spectrum(intrng(1):intrng(2),GH(3,i)));
            integral(j,i,:) = [I,GH(2,i),GH(1,i)];
        end
    end
endfunction



