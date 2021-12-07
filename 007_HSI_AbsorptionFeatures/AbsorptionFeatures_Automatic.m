function [C, Cindx] = AbsorptionFeatures_Automatic(M,RGB,d,wl)
% Mapping the wavelength position of deepest absorption features 
% to explore mineral diversity in hyperspectral images for specific
% wavelength ranges :
%   VNIR: 'Mineral','Iron','?','Iron oxydes'
%   SWIR: 'Mineral','Al-OH','Fe-OH','CO3/Mg-OH','Clay-Carbonates','Calcite-Illite'

% RGB or pseudo_RGB creation
if nargin<3
    if size(M,3)==98
        RGB=M(:,:,[50 26 13]);
    else
        RGB=M(:,:,[107 111 124]);
    end
end

if median(wl)>1000
    answer = questdlg('Would you like to reduce between 2000-2500 nm?', ...
        'Data reduction','Yes','No','Yes');
    if strcmp(answer,'Yes')
        [~,idx]=find(wl>2000);
        wl=wl(idx);
        M=M(:,:,idx);
    end
end

% Continuum removal preprocessing
tic
[Mcr,~,wl]=continuum_removal(wl,M);
toc 

% Find main absorptions
a=[];b=[];pks=[];locs=[];
h=waitbar(0,'Find main absorptions');
iter=1;
for i=1:size(Mcr,1)
    for j=1:size(Mcr,2)
        waitbar(iter/(size(Mcr,1)*size(Mcr,2)))
        locsi=[];pksi=[];
        % Smoothing with savitsky golay
        S=savgol(squeeze(Mcr(i,j,:))',7,2,0)*-1;
        % Peaks determination
        %         [pksi,locsi] = findpeaks(S);
        Sp=zeros(1,length(S));
        for l=4:size(S,2)-3
            if S(1,l)>S(1,l-3)&&S(1,l)>S(1,l+3)&&S(1,l)>S(1,l-2)&&S(1,l)>S(1,l+2)&&S(1,l)>S(1,l-1)&&S(1,l)>S(1,l+1)
                Sp(l)=1;
            else
                Sp(l)=0;
            end
        end
        [~,locsi]=find(Sp==1);
        pksi=S(1,Sp==1);
        for k=1:length(locsi)
            p=[];
            if locsi(k)>3&&locsi(k)<length(wl)-3
                % Modeling of the peak with hyperbola
                p=polyfit(1:7,S(locsi(k)-3:locsi(k)+3),2);
                % Estimation of the location of the maximum of absorption
                locs{i,j}(1,k)=locsi(k)-4-p(2)/(2*(p(1)));
                % Estimation of the maximum of absorption
                pks{i,j}(1,k)=polyval(p,-p(2)/(2*(p(1))));
            else
                locs{i,j}(1,k)=locsi(k);
                pks{i,j}(1,k)=pksi(k);
            end
        end
        if isempty(k)==0
            a=[a pks{i,j}];
            b=[b locs{i,j}];
        end
        iter=iter+1;
    end
end
close(h)

% Main absorption of all the pixels
Y = discretize(b,0:1:size(Mcr,3));

% Graphical visualization of the main absorption
figure;
histogram(Y,0:1:size(Mcr,3));
grid on;
set(gca,'Xtick',0:15:size(Mcr,3),'Xticklabel',round(wl(1:15:end)),'fontsize',14)
xlabel('Wavelength (nm)')
ylabel('Proportion')
title('Absorption features')

% Selection of the absorption areas
if median(wl)<1000
    % Position to study
    wlabs=[415 500;675 800;870 960];
    % Proxy studied
    label={'Mineral','Iron','?','Iron oxydes'};
    for i=1:size(wlabs,1)
        for j=1:size(wlabs,2)
            [a,~]=find(abs(wl-wlabs(i,j))==min(abs(wl-wlabs(i,j))));
            Cindx(i,j)=a(1);
        end
    end
else
    % Position to study
    if strcmp(answer,'Yes')
        wlabs=[2160 2220;2230 2295;2300 2360; 2160 2360;2400 2480];
        % Proxy studied
        label={'Al-OH','Fe-OH','CO3/Mg-OH','Clay-Carbonates','Calcite-Illite'};
    else
        wlabs=[1000 1170; 2160 2220;2230 2295;2300 2360; 2160 2360;2400 2480];
        % Proxy studied
        label={'Mineral','Al-OH','Fe-OH','CO3/Mg-OH','Clay-Carbonates','Calcite-Illite'};
    end
    for i=1:size(wlabs,1)
        for j=1:size(wlabs,2)
            [a,~]=find(abs(wl-wlabs(i,j))==min(abs(wl-wlabs(i,j))));
            Cindx(i,j)=a(1);
        end
    end
end

% Estimation of the mineral map (absorption and wavelength)
h=waitbar(0,'Mineral map estimations');
iter=1;
for k=1:size(Cindx,1)
    c=zeros(size(Mcr,1),size(Mcr,2),3);
    for i=1:size(Mcr,1)
        for j=1:size(Mcr,2)
            b=[];
            waitbar(iter/(size(Mcr,1)*size(Mcr,2)*size(Cindx,1)))
            if sum(locs{i,j}>Cindx(k,1)&locs{i,j}<Cindx(k,2))>0
                if sum(locs{i,j}>Cindx(k,1)&locs{i,j}<Cindx(k,2))==1
                    [~,b]=find(locs{i,j}>Cindx(k,1)&locs{i,j}<Cindx(k,2));
                    c(i,j,1)=pks{i,j}(b)*-1;
                    c(i,j,2)=1;
                    c(i,j,3)=((locs{i,j}(b)-Cindx(k,1))/(Cindx(k,2)-Cindx(k,1)))*(wl(Cindx(k,2))-wl(Cindx(k,1)))+wl(Cindx(k,1));
                else
                    [~,b]=find(locs{i,j}>Cindx(k,1)&locs{i,j}<Cindx(k,2));
                    [~,b]=find(pks{i,j}==max(pks{i,j}(b)));
                    c(i,j,1)=pks{i,j}(b(1))*-1;
                    c(i,j,2)=1;
                    c(i,j,3)=((locs{i,j}(b(1))-Cindx(k,1))/(Cindx(k,2)-Cindx(k,1)))*(wl(Cindx(k,2))-wl(Cindx(k,1)))+wl(Cindx(k,1));
                end
            else
                c(i,j,1:3)=NaN;
            end
            iter=iter+1;
        end
    end
    C{k}=c;
    
    % Display the map
    figure
    ha(1)=subplot(311); % RGB image
    imagesc(d,d(1:size(RGB,1)),RGB*(0.5/mean(RGB(:))))
    if Cindx(k,1)<1
        title(strcat(num2str(round(wl(1))),'-',num2str(round(wl(Cindx(k,2)))),' nm:',label{k}))
    else if Cindx(k,2)>length(wl)
            title(strcat(num2str(round(wl(Cindx(k,1)))),'-',num2str(round(wl(end))),' nm:',label{k}))
        else
            title(strcat(num2str(round(wl(Cindx(k,1)))),'-',num2str(round(wl(Cindx(k,2)))),' nm:',label{k}))
        end
    end
    colorbar
    ylabel('Width (cm)')
    set(gca,'fontsize',14)
    ha(2)=subplot(312); % Absorption intensity map
    imagesc(d,d(1:size(RGB,1)),1-squeeze(c(:,:,1)))
    colorbar
    caxis(sort(1-[nanmax([nanmin(reshape(squeeze(c(:,:,1)),[],1)) nanmedian(reshape(squeeze(c(:,:,1)),[],1))-3*nanstd(reshape(squeeze(c(:,:,1)),[],1))]) nanmin([nanmax(reshape(squeeze(c(:,:,1)),[],1)) nanmedian(reshape(squeeze(c(:,:,1)),[],1))+3*nanstd(reshape(squeeze(c(:,:,1)),[],1))])]))
    colormap([0 0 0; jet])
    title('Absorption intensity (normalized reflectance)')
    ylabel('Width (cm)')
    set(gca,'fontsize',14)
    ha(3)=subplot(313); % Absorption wavelength map
    imagesc(d,d(1:size(RGB,1)),squeeze(c(:,:,3)))
    colormap([1 1 1; jet])
    colorbar
    caxis([nanmax([nanmin(reshape(squeeze(c(:,:,3)),[],1)) nanmedian(reshape(squeeze(c(:,:,3)),[],1))-3*nanstd(reshape(squeeze(c(:,:,3)),[],1))]) nanmin([nanmax(reshape(squeeze(c(:,:,3)),[],1)) nanmedian(reshape(squeeze(c(:,:,3)),[],1))+3*nanstd(reshape(squeeze(c(:,:,3)),[],1))])])
    title('Absorption wavelength (relative position in the interval)')
    ylabel('Width (cm)')
    xlabel('Depth (cm)')
    set(gca,'fontsize',14)
    linkaxes(ha,'xy')
end
close(h)

end