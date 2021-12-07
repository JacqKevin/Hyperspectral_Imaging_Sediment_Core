function [C,label,clrmap] = IndiceCalculation(M,RGB,d,wl,indx,norm,fig)
% Function to calculate spectral indices
% INPUT :
%           M : Hyperspectral datacube
%           RGB : RGB image of the sample
%           d: Depth vector
%           wl : Wavelength
%           optional:
%               indx : Index of the index to estimate (0 to choose after)
%               norm : Choose if you don't want to normalize (0), or by
%                   Rmean (1) or by d555 (2, only for VNIR data)
%               fig : Diplay or not figures (1/0)
% OUTPUT :
%           C : Index map(s)
%           label : Name of the estimated indices
%           clrmap : Colormaps used for the indices

if length(wl)==98||length(wl)==144
    wl=wl(10:end-10);
    M=M(:,:,10:end-10);
end

RGB=RGB*(0.5/median(RGB(:)));

if nargin<7
    fig=1;
end

if mean(M(:))>100
    M=M/10000*100; % Reflectance
end

Rmean=mean(M,3);

bandwidth=median(wl(2:end)-wl(1:end-1));

map_red=cat(2,[0.5:0.001:1;zeros(1,501);zeros(1,501)],[ones(1,1001);0:0.001:1;0:0.001:1])';
map_green=cat(2,[zeros(1,1001);0.5:0.0005:1;zeros(1,1001)],[0:0.001:1;ones(1,1001);0:0.001:1])';
map_orange=[ones(1,1501);1:-0.0005:.25;zeros(1,1501)]';
map_blue=cat(2,[zeros(1,501);zeros(1,501);0.5:0.001:1],[0:0.001:1;0:0.001:1;ones(1,1001)])';
map_gray=[0:0.001:1;0:0.001:1;0:0.001:1]';
map_purple=[1:-0.00005:0.5;1:-0.0001:0;1:-0.00005:0.5]';
p_r=[1 0 0];
p_g=[0 1 0];
p_b=[0 0 1];
p_o=[1 1 0];
p_gr=[0.5 0.5 0.5];
p_p=[0.5 0 0.5];

%% VNIR
% Check the wavelengths range
iter=0;
if sum(wl<800&wl>600)>0
    
    index_vnir={'Chlorophylls R675/R750',...
        'Chlorophylls R645/R675',...
        'Chlorophylls R660/R670',...
        'Chlorophylls R590/R690',...
        'Chlorophylls R675-R750',...
        'Chlorophylls R645-R675',...
        'Chlorophylls dR660-dR690',...
        'Chlorophylls d675',...
        'Chlorophylls d690',...
        'Chlorophylls RABD670',...
        'Chlorophylls I-band',...
        'Chlorophylls Area650-750 (eq Wolfe)',...
        'Chlorophylls Area650-700 (eq Wolfe)',...
        'Chlorophylls Area590-730 (eq Wolfe)',...
        'Chlorophylls Area600-760 (eq Wolfe)',...
        'Chlorophylls Area650-750 (eq Trachsel)',...
        'Chlorophylls Area650-700 (eq Trachsel)',...
        'Chlorophylls Area590-730 (eq Trachsel)',...
        'Chlorophylls Area600-760 (eq Trachsel)',...
        'Chlorophylls (Area650-700)/R670',...
        'Phycocyanin (Area600-630)/R615',...
        'Bacteriochlorophyll a RABD846',...
        'Bacteriochlorophyll a Area750-900/R844',...
        'Carotenoids RABD510',...
        'Carotenoids RABD550',...
        'Total pigment content RABA400-560',...
        'Oxydes d555',...
        'Oxydes d575 (Hematite)',...
        'Oxydes d535 (Goethite)',...
        'Hm/Gt d555/d535',...
        'Iron Oxydes R720/R880'};

    associated_vnir_map={map_green;...
        map_green;...
        flip(map_green);...
        map_green;...
        map_green;...
        map_green;...
        map_green;...
        flip(map_green);...
        map_green;...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_green);...
        flip(map_blue);...
        map_purple;...
        map_purple;...
        map_orange;...
        map_orange;...
        jet(1500);...
        flip(map_red);...
        map_red;...
        map_red;...
        map_red;...
        flip(map_red)};

    associated_vnir_p={p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        p_g;...
        [0 0 1];...
        p_p;...
        p_p;...
        p_o;...
        p_o;...
        [0 0 0];...
        p_r;...
        p_r;...
        p_r;...
        p_r;...
        p_r};

    if nargin<5||(length(indx)==1&&indx==0)
        [indx,~] = listdlg('ListString',index_vnir,'ListSize',[250 300]);
    end
    
    for i=1:length(indx)
        index_vnir_selec{i}=index_vnir{indx(i)};
    end
    
    if sum(double(strcmp('Chlorophylls R675/R750',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[675 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls R645/R675',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[645 675];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls R660/R670',index_vnir_selec)))
        %von Gunten, L., Grosjean, M. (2009) High-resolution, quantitative climate reconstruction over the past 1000 years and pollution history derived from lake sediments in Central Chile. Philosophisch-naturwissenschaftlichen Fakultät 246
        iter=iter+1;
        wlr=[660 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls R590/R690',index_vnir_selec)))
        %Trachsel, M., Grosjean, M., Schnyder, D., Kamenik, C., Rein, B. (2010) Scanning reflectance spectroscopy (380–730 nm): a novel method for quantitative high-resolution climate reconstructions from minerogenic lake sediments. Journal of Paleolimnology 44: 979–994
        iter=iter+1;
        wlr=[590 690];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls R675-R750',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[675 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls R645-R675',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[645 675];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls dR660-dR690',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[660 690];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        
        IM=squeeze(Md(:,:,ind1(1)))-squeeze(Md(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Chlorophylls d675',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=675;
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Chlorophylls d690',index_vnir_selec)))
        %Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=690;
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Chlorophylls RABD670',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[590 730 660 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        [~,ind4]=find(abs(wl-wlr(1,4))==min(abs(wl-wlr(1,4))));
        
        [t,I] = min(M(:,:,ind3(1):ind4(1)),[],3);
        for i=1:size(M,1)
            for j=1:size(M,2)
                IM(i,j)=(((ind2(1)-ind3(1)+I(i,j)-1)*squeeze(M(i,j,ind1(1)))+...
                    (ind3(1)+I(i,j)-1-ind1(1))*squeeze(M(i,j,ind2(1))))/...
                    (ind2(1)-ind1(1)))./t(i,j);
            end
        end
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls I-band',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[590 730 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((ind2(1)-ind3(1))*squeeze(M(:,:,ind1(1)))+...
            (ind3(1)-ind1(1))*squeeze(M(:,:,ind2(1))))/...
            (ind2(1)-ind1(1)))./squeeze(M(:,:,ind3(1)))./Rmean;
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area650-750 (eq Wolfe)',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area650-700 (eq Wolfe)',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area590-730 (eq Wolfe)',index_vnir_selec)))
        % Trachsel, M., Grosjean, M., Schnyder, D., Kamenik, C., Rein, B. (2010) Scanning reflectance spectroscopy (380–730 nm): a novel method for quantitative high-resolution climate reconstructions from minerogenic lake sediments. Journal of Paleolimnology 44: 979–994
        iter=iter+1;
        wlr=[590 730];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area600-760 (eq Wolfe)',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[600 760];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area650-750 (eq Trachsel)',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area650-700 (eq Trachsel)',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area590-730 (eq Trachsel)',index_vnir_selec)))
        % Trachsel, M., Grosjean, M., Schnyder, D., Kamenik, C., Rein, B. (2010) Scanning reflectance spectroscopy (380–730 nm): a novel method for quantitative high-resolution climate reconstructions from minerogenic lake sediments. Journal of Paleolimnology 44: 979–994
        iter=iter+1;
        wlr=[590 730];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls Area600-760 (eq Trachsel)',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[600 760];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Chlorophylls (Area650-700)/R670',index_vnir_selec)))
        % Van Exem
        iter=iter+1;
        wlr=[650 700 670];
        %         wlr=[630 700 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        %                 IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
        %                     (ind3-ind1)*squeeze(M(:,:,ind2(1))))/...
        %                     (ind2-ind1))./squeeze(M(:,:,ind3(1))); %RABD
        %                 IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
        %                     ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
        %                     bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
        %                     bandwidth/2*squeeze(M(:,:,ind2(1)))-...
        %                     bandwidth/2*squeeze(M(:,:,ind1(1))))./...
        %                     (squeeze(M(:,:,ind3(1))));% RABA/Area Trachsel
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            squeeze(sum(M(:,:,(ind1(1)):(ind2(1))),3)))./...
            (squeeze(M(:,:,ind3(1))));% RABA/Area Das
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Phycocyanin (Area600-630)/R615',index_vnir_selec)))
        % Sorrel, P., Jacq, K., et al. (2021) Evidence for centennial-scale Mid-Holocene episodes of hypolimnetic anoxia in a high-altitude lake system from central Tian Shan (Kyrgyzstan). Quaternary Science Reviews 252:
        iter=iter+1;
        wlr=[600 630 615];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))))./...
            (squeeze(M(:,:,ind3(1))));% RABA/Area
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Bacteriochlorophyll a RABD846',index_vnir_selec)))
        % Butz, C., Grosjean, M., Fischer, D., Wunderle, S., Tylmann, W., Rein, B. (2015) Hyperspectral imaging spectroscopy: a promising method for the biogeochemical analysis of lake sediments. Journal of Applied Remote Sensing 9: 1–20
        iter=iter+1;
        wlr=[790 900 845];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        %                 IM=((34*squeeze(M(:,:,ind1(1)))+...
        %                     34*squeeze(M(:,:,ind2(1))))/...
        %                     (34+34))./squeeze(M(:,:,ind3(1)));
        IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
            (ind3-ind1)*squeeze(M(:,:,ind2(1))))/...
            (ind2-ind1))./squeeze(M(:,:,ind3(1))); %RABD
        
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Bacteriochlorophyll a Area750-900/R844',index_vnir_selec)))
        % Sorrel, P., Jacq, K., et al. (2021) Evidence for centennial-scale Mid-Holocene episodes of hypolimnetic anoxia in a high-altitude lake system from central Tian Shan (Kyrgyzstan). Quaternary Science Reviews 252:
        iter=iter+1;
        wlr=[750 900 844];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))))./...
            (squeeze(M(:,:,ind3(1))));% RABA/Area
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Carotenoids RABD510',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[490 530 510];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
            (ind3-ind1)*squeeze(M(:,:,ind2(1))))/2)./squeeze(M(:,:,ind3(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Carotenoids RABD550',index_vnir_selec))) % A supprimer?
        % Van Exem
        iter=iter+1;
        wlr=[538 576 550];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        %         IM=((34*squeeze(M(:,:,ind1(1)))+...
        %             34*squeeze(M(:,:,ind2(1))))/...
        %             (34+34))./squeeze(M(:,:,ind3(1)));
        IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
            (ind3-ind1)*squeeze(M(:,:,ind2(1))))/...
            (ind2-ind1))./squeeze(M(:,:,ind3(1))); %RABD
        %         IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
        %             ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2-...
        %             bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
        %             bandwidth/2*squeeze(M(:,:,ind2(1)))-...
        %             bandwidth/2*squeeze(M(:,:,ind1(1))))./...
        %             (squeeze(M(:,:,ind3(1))));% RABA/Area
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Total pigment content RABA400-560',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[400 560 590];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sum(repmat(squeeze(M(:,:,ind3(1))),1,1,ind2(1)-ind1(1)+1)./M(:,:,ind1(1):ind2(1)),3);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Oxydes d555',index_vnir_selec)))
        % Deaton, B.C., Balsam, W.L., 1991. Visible spectroscopy - a rapid method for determining hematite and goethite concentration in geological materials. J. Sediment. Petrol. 61, 628–632. https://doi.org/10.1306/d4267794-2b26-11d7-8648000102c1865d
        iter=iter+1;
        wlr=555;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Oxydes d575 (Hematite)',index_vnir_selec)))
        % Deaton, B.C., Balsam, W.L., 1991. Visible spectroscopy - a rapid method for determining hematite and goethite concentration in geological materials. J. Sediment. Petrol. 61, 628–632. https://doi.org/10.1306/d4267794-2b26-11d7-8648000102c1865d
        iter=iter+1;
        wlr=575;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Oxydes d535 (Goethite)',index_vnir_selec)))
        % Deaton, B.C., Balsam, W.L., 1991. Visible spectroscopy - a rapid method for determining hematite and goethite concentration in geological materials. J. Sediment. Petrol. 61, 628–632. https://doi.org/10.1306/d4267794-2b26-11d7-8648000102c1865d
        iter=iter+1;
        wlr=535;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Hm/Gt d555/d535',index_vnir_selec)))
        %Wu, G., Xu, T., Zhang, X., Zhang, C., Yan, N., 2016. The visible spectroscopy of iron oxide minerals in dust particles from ice cores on the Tibetan Plateau. Tellus, Ser. B Chem. Phys. Meteorol. 68, 1–10. https://doi.org/10.3402/tellusb.v68.29191
        iter=iter+1;
        wlr=[555 535];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7*(length(wl)/98),2,1),size(M,1),size(M,2),size(M,3));
        IM=squeeze(Md(:,:,ind1(1)))./squeeze(Md(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
        
        clear Md
    end
    
    if sum(double(strcmp('Iron Oxydes R720/R880',index_vnir_selec)))
        % Jackisch, R., Lorenz, S., Zimmermann, R., Möckel, R., Gloaguen, R. (2018) Drone-borne hyperspectral monitoring of acid mine drainage: An example from the Sokolov lignite district. Remote Sensing 10: 1–23
        iter=iter+1;
        wlr=[720 880];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if nargin<6
        Norm = questdlg('Would you like to normalize the indices?', ...
            'Normalization', 'Rmean','d555','No','No');
    else
        if norm==0
            Norm='No';
        else
            if norm==1
                Norm='Rmean';
            else
                if norm==2
                    Norm='d555';
                end
            end
        end
    end
    
    if strcmp(Norm,'Rmean')
        for i=1:size(C,3)
            IM=C(:,:,i)./Rmean;
            [a,b]=find(IM==inf);
            for j=1:length(a)
                IM(a(j),b(j))=0;
            end
            [a,b]=find(IM==-inf);
            for j=1:length(a)
                IM(a(j),b(j))=0;
            end
            [a,b]=find(isnan(IM)==1);
            for j=1:length(a)
                if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                    IM(a(j),b(j))=nanmedian(IM(:));
                else
                    IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
                end
            end
            C(:,:,i)=IM;
        end
    else if strcmp(Norm,'d555')
            D=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
            [~,ind]=find(abs(wl-555)==min(abs(wl-555)));
            d555=squeeze(D(:,:,ind));
            for i=1:size(C,3)
                IM=squeeze(C(:,:,i))./d555;
                [a,b]=find(IM==inf);
                for j=1:length(a)
                    IM(a(j),b(j))=0;
                end
                [a,b]=find(IM==-inf);
                for j=1:length(a)
                    IM(a(j),b(j))=0;
                end
                [a,b]=find(isnan(IM)==1);
                for j=1:length(a)
                    if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                        IM(a(j),b(j))=nanmedian(IM(:));
                    else
                        IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
                    end
                end
                C(:,:,i)=IM;
            end
        end
    end
    
    if fig==1
        if sum(wl<2400&wl>1100)==0
            
            NbSubPerFig=4; % 5 subplot
            
            NumFig=ceil(size(C,3)/NbSubPerFig);
            for j=1:NumFig
                figure
                ha(1)=subplot(1,NbSubPerFig+1,1);
                imagesc(d(1:size(IM,1))-d(1),d,imrotate(RGB,-90))
                ylabel('Depth (cm)')
                colorbar('southoutside');
                set(ha(1),'fontsize',14)
                for i=1:NbSubPerFig
                    if i+(j-1)*NbSubPerFig<=size(C,3)
                        IM=squeeze(C(:,:,i+(j-1)*NbSubPerFig));
                        ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                        imagesc(d(1:size(IM,1))-d(1),d,imrotate(squeeze(C(:,:,i+(j-1)*NbSubPerFig)),-90))
                        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
                        colorbar('southoutside')
                        map=associated_vnir_map{indx(i+(j-1)*NbSubPerFig)};
                        colormap(ha(i+1),map)
                        clrmap{i+(j-1)*NbSubPerFig}=map;
                        %         xlabel('Width (cm)')
                        if i==NbSubPerFig
                            ylabel('Depth (cm)')
                        end
                        title(index_vnir{indx(i+(j-1)*NbSubPerFig)})
                        if i==NbSubPerFig
                            set(ha(i+1),'fontsize',14,'YAxisLocation','right');
                        else
                            set(ha(i+1),'fontsize',14,'ytick',[]);
                        end
                    else
                        ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                        imagesc(d(1:size(IM,1))-d(1),d,imrotate(zeros(size(C,1),size(C,2)),-90))
                        colorbar('southoutside')
                        set(ha(i+1),'fontsize',14,'ytick',[]);
                    end
                end
                linkaxes(ha,'xy')
                
                figure
                ha(1)=subplot(1,NbSubPerFig+1,1);
                imagesc(d(1:size(IM,1))-d(1),d,imrotate(RGB,-90))
                ylabel('Depth (cm)')
                set(ha(1),'fontsize',14)
                for i=1:NbSubPerFig
                    if i+(j-1)*NbSubPerFig<=size(C,3)
                        IM=squeeze(C(:,:,i+(j-1)*NbSubPerFig));
                        ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                        plot(nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),d,'color',associated_vnir_p{indx(i+(j-1)*NbSubPerFig)})
                        grid on
                        grid minor
                        ylim([min(d) max(d)])
                        if i==NbSubPerFig
                            ylabel('Depth (cm)')
                        end
                        title(index_vnir{indx(i+(j-1)*NbSubPerFig)})
                        if i==NbSubPerFig
                            set(ha(i+1),'fontsize',14,'YAxisLocation','right','ydir','reverse');
                        else
                            set(ha(i+1),'fontsize',14,'ytick',[],'ydir','reverse');
                        end
                    end
                end
                linkaxes(ha,'y')
            end
            
            figure;
            imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
            colorbar
            colormap(jet)
            xtickangle(45)
            set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
        else
            labelvnir=label;
        end
    else
        for i=1:size(C,3)
            clrmap{i}=associated_vnir_map{indx(i)};
        end
    end
end

if sum(wl<2400&wl>1100)>0
    
    index_swir={'OM R1650',...
        'Aromatic Organic Matter',...
        'Clay R2200',...
        'Clay RABA2190-2230',...
        'Hydroxyl bonds R1450',...
        'Moisture R1935',...
        'Moisture R1450/R1300',...
        'R1935/R1450',...
        'Normalized Differenced Gypsum Ratio (NDGI)',...
        'Calcite R2340',...
        'Illite-muscovite Crystallinity R1900/R2200',...
        'Kaolinite Crystallinity Index (KCI)',...
        'TBC IND2274',...
        'TBC IND1754',...
        'TBC Ratio2396',...
        'TBC MLR1',...
        'TBC MLR2',...
        'Senna Index (Kaolinite Crystallinity)',...
        '22SP-Index 2200nm',...
        'n smectites alumineuses',...
        'n gypse'};

    associated_swir_map={flip(jet(1500));...
        jet(1500);...
        map_orange;...
        map_orange;...
        flip(map_blue);...
        flip(map_blue);...
        flip(map_blue);...
        flip(map_blue);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        jet(1500);...
        map_orange;...
        jet(1500);...
        jet(1500)};
    
    associated_swir_p={[0 0 0];...
        [0 0 0];...
        p_o;...
        p_o;...
        p_b;...
        p_b;...
        p_b;...
        p_b;...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        [0 0 0];...
        p_o;...
        [0 0 0];...
        [0 0 0]};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_swir);
    end
    
    for i=1:length(indx)
        index_swir_selec{i}=index_swir{indx(i)};
    end
    
    if sum(double(strcmp('OM R1650',index_swir_selec)))
        iter=iter+1;
        wlr=1650;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Aromatic Organic Matter',index_swir_selec)))
        iter=iter+1;
        wlr=[1660 1690 1670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        %         IM=((34*squeeze(M(:,:,ind1(1)))+...
        %             34*squeeze(M(:,:,ind2(1))))/...
        %             (34+34))./squeeze(M(:,:,ind3(1)));
        %         IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
        %             (ind3-ind1)*squeeze(M(:,:,ind2(1))))/...
        %             (ind2-ind1))./squeeze(M(:,:,ind3(1))); %RABD
        IM=(squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))+...
            ((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*...
            (wl(ind2(1))-wl(ind1(1))))/2-...
            bandwidth*squeeze(sum(M(:,:,(ind1(1)+1):(ind2(1)-1)),3))-...
            bandwidth/2*squeeze(M(:,:,ind2(1)))-...
            bandwidth/2*squeeze(M(:,:,ind1(1))))./...
            (squeeze(M(:,:,ind3(1))));% RABA/Area
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Clay R2200',index_swir_selec)))
        iter=iter+1;
        wlr=2200;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Clay RABA2190-2230',index_swir_selec)))
        iter=iter+1;
        wlr=[2190 2230];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)))*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Hydroxyl bonds R1450',index_swir_selec)))
        iter=iter+1;
        wlr=1450;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Moisture R1935',index_swir_selec)))
        iter=iter+1;
        wlr=1935;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Moisture R1450/R1300',index_swir_selec)))
        % Bull, C. R. (1991). Wavelength selection for near-infrared reflectance moisture meters. Journal of Agricultural Engineering Research, 49(C), 113–125. https://doi.org/10.1016/0021-8634(91)80032-A
        iter=iter+1;
        wlr=[1450 1300];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1))))./(squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('R1935/R1450',index_swir_selec)))
        iter=iter+1;
        wlr=[1935 1450];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1))))./(squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Differenced Gypsum Ratio (NDGI)',index_swir_selec)))
        % Milewski, R., Chabrillat, S., Brell, M., Schleicher, A.M., Guanter, L. (2019) Assessment of the 1.75 µm absorption feature for gypsum estimation using laboratory, air- and spaceborne hyperspectral sensors. International Journal of Applied Earth Observation and Geoinformation 77: 69–83
        iter=iter+1;
        wlr=[1690 1750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('n gypse',index_swir_selec)))
        % Stage Marqueurs spectroscopiques des changements paléoenvironnementaux en Mer d’Aral à l’Holocène récent
        iter=iter+1;
        wlr=[1750 1710 1780];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        IM=(1-squeeze(M(:,:,ind1(1))))./...
            (squeeze(M(:,:,ind2(1)))*((wlr(3)-wlr(1))/(wlr(3)-wlr(2)))+...
            squeeze(M(:,:,ind3(1)))*((wlr(1)-wlr(2))/(wlr(3)-wlr(2))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('n smectites alumineuses',index_swir_selec)))
        % Stage Marqueurs spectroscopiques des changements paléoenvironnementaux en Mer d’Aral à l’Holocène récent
        iter=iter+1;
        wlr=[2210 2179 2248];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        IM=(1-squeeze(M(:,:,ind1(1))))./...
            (squeeze(M(:,:,ind2(1)))*((wlr(3)-wlr(1))/(wlr(3)-wlr(2)))+...
            squeeze(M(:,:,ind3(1)))*((wlr(1)-wlr(2))/(wlr(3)-wlr(2))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Calcite R2340',index_swir_selec)))
        % Sun, L., Khan, S., Godet, A. (2018) Integrated ground-based hyperspectral imaging and geochemical study of the Eagle Ford Group in West Texas. Sedimentary Geology 363: 34–47
        iter=iter+1;
        wlr=2340;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('TBC IND2274',index_swir_selec)))
        % Rivard, B., Lyder, D., Feng, J., Gallie, A., Cloutis, E., Dougan, P., Gonzalez, S., Cox, D., Lipsett, M.G., 2010. Bitumen content estimation of Athabasca oil sand from broad band infrared reflectance spectra. Can. J. Chem. Eng. 88, 830–838. https://doi.org/10.1002/cjce.20343
        iter=iter+1;
        wlr=[2210 2274];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    
    if sum(double(strcmp('TBC IND1754',index_swir_selec)))
        % Rivard, B., Lyder, D., Feng, J., Gallie, A., Cloutis, E., Dougan, P., Gonzalez, S., Cox, D., Lipsett, M.G., 2010. Bitumen content estimation of Athabasca oil sand from broad band infrared reflectance spectra. Can. J. Chem. Eng. 88, 830–838. https://doi.org/10.1002/cjce.20343
        iter=iter+1;
        wlr=[2054 1754];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('TBC Ratio2396',index_swir_selec)))
        % Rivard, B., Lyder, D., Feng, J., Gallie, A., Cloutis, E., Dougan, P., Gonzalez, S., Cox, D., Lipsett, M.G., 2010. Bitumen content estimation of Athabasca oil sand from broad band infrared reflectance spectra. Can. J. Chem. Eng. 88, 830–838. https://doi.org/10.1002/cjce.20343
        iter=iter+1;
        wlr=[2210 2396];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('TBC MLR1',index_swir_selec)))
        % Rivard, B., Lyder, D., Feng, J., Gallie, A., Cloutis, E., Dougan, P., Gonzalez, S., Cox, D., Lipsett, M.G., 2010. Bitumen content estimation of Athabasca oil sand from broad band infrared reflectance spectra. Can. J. Chem. Eng. 88, 830–838. https://doi.org/10.1002/cjce.20343
        iter=iter+1;
        wlr=[2220 2330];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=7.03+2.03*squeeze(M(:,:,ind1(1)))-...
            3.89*squeeze(M(:,:,ind2(1)))+...
            11.56*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('TBC MLR2',index_swir_selec)))
        % Rivard, B., Lyder, D., Feng, J., Gallie, A., Cloutis, E., Dougan, P., Gonzalez, S., Cox, D., Lipsett, M.G., 2010. Bitumen content estimation of Athabasca oil sand from broad band infrared reflectance spectra. Can. J. Chem. Eng. 88, 830–838. https://doi.org/10.1002/cjce.20343
        iter=iter+1;
        wlr=[2210 2396 2220 2330];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        [~,ind4]=find(abs(wl-wlr(1,4))==min(abs(wl-wlr(1,4))));
        IM=12.25-12.55*(squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1))))+...
            41.21*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind4(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Illite-muscovite Crystallinity R1900/R2200',index_swir_selec)))
        % Doublier, M.P., Roache, T., Potel, S. (2010) Short-wavelength infrared spectroscopy: A new petrological tool in low-grade to very low-grade pelites. Geology 38: 1031–1034
        iter=iter+1;
        wlr=[2200 1900];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1))))./(squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Kaolinite Crystallinity Index (KCI)',index_swir_selec)))
        % Alonso de Linaje, V., Khan, S.D. (2017) Mapping of diagenetic processes in sandstones using imaging spectroscopy: A case study of the Utrillas Formation, Burgos, Spain. Sedimentary Geology 353: 114–124
        iter=iter+1;
        wlr=[2120 2430 2159 2165];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [Mcr,~,wlcr]=continuum_removal(wl(ind1(1):ind2(1)),M(:,:,ind1(1):ind2(1)));
        wlcr=wlcr';
        [~,ind3]=find(abs(wlcr-wlr(1,3))==min(abs(wlcr-wlr(1,3))));
        [~,ind4]=find(abs(wlcr-wlr(1,4))==min(abs(wlcr-wlr(1,4))));
        if ind4==ind3
            ind4=ind4+1;
        end
        IM=squeeze(Mcr(:,:,ind3(1)))./squeeze(Mcr(:,:,ind4(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Senna Index (Kaolinite Crystallinity)',index_swir_selec)))
        % de Senna, J.A., de Souza Filho, C.R., Angélica, R.S. (2008) Characterization of clays used in the ceramic manufacturing industry by reflectance spectroscopy: An experiment in the São Simão ball-clay deposit, Brazil. Applied Clay Science 41: 85–98
        iter=iter+1;
        wlr=[2180 2162 2206];
        [Mcr,~,wlcr]=continuum_removal(wl,M);
        wlcr=wlcr';
        [~,ind1]=find(abs(wlcr-wlr(1,1))==min(abs(wlcr-wlr(1,1))));
        [~,ind2]=find(abs(wlcr-wlr(1,2))==min(abs(wlcr-wlr(1,2))));
        [~,ind3]=find(abs(wlcr-wlr(1,3))==min(abs(wlcr-wlr(1,3))));
        if ind1==ind2
            ind1=ind1-1;
        end
        IM=squeeze(Mcr(:,:,ind1(1)))./squeeze(Mcr(:,:,ind2(1)))+...
            1-squeeze(Mcr(:,:,ind3(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('22SP-Index 2200nm',index_swir_selec)))
        % Awad, M.E., Amer, R., López-Galindo, A., El-Rahmany, M.M., García del Moral, L.F., Viseras, C. (2018) Hyperspectral remote sensing for mapping and detection of Egyptian kaolin quality. Applied Clay Science 160: 249–262
        iter=iter+1;
        [Mcr,~,wlcr]=continuum_removal(wl,M);
        wlr=[2200 2165];
        wlcr=wlcr';
        [~,ind1]=find(abs(wlcr-wlr(1,1))==min(abs(wlcr-wlr(1,1))));
        [~,ind2]=find(abs(wlcr-wlr(1,2))==min(abs(wlcr-wlr(1,2))));
        
        IM=(squeeze(Mcr(:,:,ind1(1)))-squeeze(Mcr(:,:,ind2(1))))./...
            squeeze(Mcr(:,:,ind1(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if nargin<6
        Norm = questdlg('Would you like to normalize the indices?', ...
            'Normalization', 'Rmean','No','No');
    else
        if norm==0
            Norm='No';
        else
            if norm==1
                Norm='Rmean';
            end
        end
    end
    
    if strcmp(Norm,'Rmean')
        for i=1:size(C,3)
            IM=C(:,:,i)./Rmean;
            [a,b]=find(IM==inf);
            for j=1:length(a)
                IM(a(j),b(j))=0;
            end
            [a,b]=find(IM==-inf);
            for j=1:length(a)
                IM(a(j),b(j))=0;
            end
            [a,b]=find(isnan(IM)==1);
            for j=1:length(a)
                if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                    IM(a(j),b(j))=nanmedian(IM(:));
                else
                    IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
                end
            end
            C(:,:,i)=IM;
        end
    end
    
    if fig==1
        if sum(wl<800&wl>600)>0
            labeltmp=label;
            label=labelvnir;
            for i=1:length(labeltmp)
                label{i+length(labelvnir)}=labeltmp{i};
            end
        end
        
        NbSubPerFig=4; % 5 subplot
        
        NumFig=ceil(size(C,3)/NbSubPerFig);
        for j=1:NumFig
            figure
            ha(1)=subplot(1,NbSubPerFig+1,1);
            imagesc(d(1:size(IM,1))-d(1),d,imrotate(RGB,-90))
            ylabel('Depth (cm)')
            colorbar('southoutside');
            set(ha(1),'fontsize',14)
            for i=1:NbSubPerFig
                if i+(j-1)*NbSubPerFig<=size(C,3)
                    IM=squeeze(C(:,:,i+(j-1)*NbSubPerFig));
                    ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                    imagesc(d(1:size(IM,1))-d(1),d,imrotate(squeeze(C(:,:,i+(j-1)*NbSubPerFig)),-90))
                    caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
                    colorbar('southoutside')
                    map=associated_swir_map{indx(i+(j-1)*NbSubPerFig)};
                    colormap(ha(i+1),map)
                    clrmap{i+(j-1)*NbSubPerFig}=map;
                    %         xlabel('Width (cm)')
                    if i==NbSubPerFig
                        ylabel('Depth (cm)')
                    end
                    title(index_swir{indx(i+(j-1)*NbSubPerFig)})
                    if i==NbSubPerFig
                        set(ha(i+1),'fontsize',14,'YAxisLocation','right');
                    else
                        set(ha(i+1),'fontsize',14,'ytick',[]);
                    end
                else
                    ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                    imagesc(d(1:size(IM,1))-d(1),d,imrotate(zeros(size(C,1),size(C,2)),-90))
                    colorbar('southoutside')
                    set(ha(i+1),'fontsize',14,'ytick',[]);
                end
            end
            linkaxes(ha,'xy')
            
            figure
            ha(1)=subplot(1,NbSubPerFig+1,1);
            imagesc(d(1:size(IM,1))-d(1),d,imrotate(RGB,-90))
            ylabel('Depth (cm)')
            set(ha(1),'fontsize',14)
            for i=1:NbSubPerFig
                if i+(j-1)*NbSubPerFig<=size(C,3)
                    IM=squeeze(C(:,:,i+(j-1)*NbSubPerFig));
                    ha(i+1)=subplot(1,NbSubPerFig+1,i+1);
                    plot(nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),d,'color',associated_swir_p{indx(i+(j-1)*NbSubPerFig)})
                    grid on
                    grid minor
                    ylim([min(d) max(d)])
                    if i==NbSubPerFig
                        ylabel('Depth (cm)')
                    end
                    title(index_swir{indx(i+(j-1)*NbSubPerFig)})
                    if i==NbSubPerFig
                        set(ha(i+1),'fontsize',14,'YAxisLocation','right','ydir','reverse');
                    else
                        set(ha(i+1),'fontsize',14,'ytick',[],'ydir','reverse');
                    end
                end
            end
            linkaxes(ha,'y')
        end
        
        figure;
        imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
        colorbar
        colormap(jet)
        xtickangle(45)
        set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
    else
        for i=1:size(C,3)
            clrmap{i}=associated_swir_map{indx(i)};
        end
    end
end

if sum(wl<800&wl>600)>0&&sum(wl<2400&wl>1100)>0
    
    index_vnir_swir={'SOM'};
    associated_vnir_swir_map={map_green};
    associated_vnir_swir_p={p_g};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_vnir_swir);
    end
    
    for i=1:length(indx)
        index_vnir_swir_selec{i}=index_vnir_swir{indx(i)};
    end
    
    if sum(double(strcmp('SOM',index_vnir_swir_selec)))
        % Tian, Y., Zhang, J., Yao, X., Cao, W., Zhu, Y. (2013) Laboratory assessment of three quantitative methods for estimating the organic matter content of soils in China based on visible/near-infrared reflectance spectra. Geoderma 202–203: 161–170
        iter=iter+1;
        wlr=[554 1398];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
        
        IM=ln(squeeze(Md(:,:,ind1(1)))-squeeze(Md(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            IM(a(j),b(j))=0;
        end
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=nanmedian(IM(:));
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir_swir{indx(iter)};
    end
    
    for i=1:size(C,3)
        ha=[];
        IM=squeeze(C(:,:,i));
        figure;
        ha(1)=subplot(3,1,1);
        imagesc(d,d(1:size(C,1)),IM)
        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
        colorbar
        map=associated_vnir_swir_map{indx(i)};
        clrmap{i}=map;
        colormap(map)
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        title(index_vnir_swir{indx(i)})
        ha(2)=subplot(3,1,2);
        imagesc(d,d(1:size(IM,1)),RGB)
        colorbar
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        ha(3)=subplot(313);
        plot(d,nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),'color',associated_vnir_swir_p{indx(i)})
        grid on
        xlabel('Depth (cm)')
        ylabel('SI')
        set(gca,'fontsize',14)
        colorbar
        xlim([min(d) max(d)])
        linkaxes(ha,'x')
    end
    
    if sum(wl<800&wl>600)>0
        labeltmp=label;
        label=labelvnir;
        for i=1:length(labeltmp)
            label{i+length(labelvnir)}=labeltmp{i};
        end
    end
    
    figure;
    imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
    colorbar
    colormap(jet)
    xtickangle(45)
    set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
end
end