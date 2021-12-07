function HSIvisualize(M,RGB,d,wl,indx)
% Function to visualize the HSI cube based on :
%      * Composite images
%      * Abberant pixels
%      * Median and standard deviation spectra
%      * Raw, Continuum removed, Continuum and First derivative
%            spectral profiles along the sample
%      * Histogram
%      * Wavelength correlations
%      * Q7/4 vs L* diagram
%      * PCA
%      * Minimum Noise Fraction compression to create a composite image
%      * Clustering
%      * Texture analysis
% INPUT:
%           M: Hyperspectral datacube (n*m*p)
%           RGB: Associated RGB image (n*m*3)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
%           indx: method to use:
%               0: all
%               1: Images (RGB, CIR, NIR, L*a*b*, HSV) ou (pRGB, hypRGB, HC1, HC2, L*a*b*, HSV)
%               2: 'Aberrant pixels'
%               3: 'SNR (local)' -> voir pour flou
%               4: 'Median and standard deviation spectra'
%               5: '10 most different spectra (Kennard and Stone)'
%               6: 'Raw spectra (2D)'
%               7: 'Raw spectra (3D)'
%               8: 'Continuum removed (2D)'
%               9: 'Continuum removed (2D enhanced)'
%               10: 'Continuum removed (3D)'
%               11: 'Continuum (2D)'
%               12: 'Continuum (3D)'
%               13: 'FDS (2D)'
%               14: 'FDS (3D)'
%               15: 'Wavelength correlation'
%               16: 'Q700/500 vs L* (global)'
%               17: 'Q700/500 vs L* (central)'
%               18: 'Q700/500 vs L* (map)'
%               19: 'Q700/500 vs L* with classification'
%               20: 'Principal Components Analysis (PCA)'
%               21: 'Maximum Noise Fraction (MNF) + Superpixels'
%               22: 'Clustering (Kmeans, HAC)'
%               23: 'Texture analysis'

%If M is in reflectance *10000
m=mean(M(:));
if m>100
    M=M/10000*100;
end

% RGB enhancement
RGB=RGB*(0.5/mean(RGB(:)));

% Check unique depth
[~,ia,~] = unique(d);

if median(wl)<1000
    list = {'Images (RGB, CIR, NIR, L*a*b*, HSV)',...
        'Aberrant pixels',...
        'SNR (local)',...
        'Median and standard deviation spectra',...
        '10 most different spectra (Kennard and Stone)',...
        'Grayscale histogram',...
        'Grayscale histogram as a function of depth',...
        'Spectral histogram',...
        'Spectral histogram as a function of depth',...
        'Raw spectra (2D)',...
        'Raw spectra (3D)',...
        'Continuum removed (2D)',...
        'Continuum removed (2D enhanced)',...
        'Continuum removed (3D)',...
        'Continuum (2D)',...
        'Continuum (3D)',...
        'FDS (2D)',...
        'FDS (3D)',...
        'Wavelength correlation',...
        'Q700/500 vs L* (global)',...
        'Q700/500 vs L* (central)',...
        'Q700/500 vs L* (map)',...
        'Q700/500 vs L* with classification',...
        'Principal Component Analysis (PCA)',...
        'Minimum Noise Fraction (MNF) + Superpixels',...
        'Clustering (Kmeans, HAC)',...
        'Texture analysis'};
else
    list = {'Images (pRGB, hypRGB, HC1, HC2, L*a*b*, HSV)',...
        'Aberrant pixels',...
        'SNR (local)',...
        'Median and standard deviation spectra',...
        '10 most different spectra (Kennard and Stone)',...
        'Grayscale histogram',...
        'Grayscale histogram as a function of depth',...
        'Spectral histogram',...
        'Spectral histogram as a function of depth',...
        'Raw spectra (2D)',...
        'Raw spectra (3D)',...
        'Continuum removed (2D)',...
        'Continuum removed (2D enhanced)',...
        'Continuum removed (3D)',...
        'Continuum (2D)',...
        'Continuum (3D)',...
        'FDS (2D)',...
        'FDS (3D)',...
        'Wavelength correlation',...
        'Principal Component Analysis (PCA)',...
        'Minimum Noise Fraction (MNF) + Superpixels',...
        'Clustering (Kmeans, HAC)',...
        'Texture analysis'};
end

if nargin<5
    [indx,~] = listdlg('ListString',list);
else
    if indx==0
        if median(wl)<1000
            indx=1:23;
        else
            indx=1:19;
        end
    end
end

for i=1:length(indx)
    index_selec{i}=list{indx(i)};
end

if sum(double(strcmp('Images (RGB, CIR, NIR, L*a*b*, HSV)',index_selec)))
    %% Images
    [RGB,R700_500,L]=Spec2rgb2Q75L(M/100,d,wl);
    assignin('base','RGB',RGB);
    assignin('base','R700_500',R700_500);
    assignin('base','L',L);
    
    figure
    ha(1)=subplot(171);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90)) % 700 550 470
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    title('RGB')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(172);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(rgb2hsv(RGB),-90))
    xlabel('Width (cm)')
    title('HSV')
    set(gca,'fontsize',14)
    
    ha(3)=subplot(173);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(rgb2lab(RGB)/100,-90))
    xlabel('Width (cm)')
    title('LAB')
    set(gca,'fontsize',14)
    
    ha(4)=subplot(174);
    % Butz, C., Grosjean, M., Fischer, D., Wunderle, S., Tylmann, W., & Rein, B. (2015). Hyperspectral imaging spectroscopy: a promising method for the biogeochemical analysis of lake sediments. Journal of Applied Remote Sensing, 9, 1–20. https://doi.org/10.1117/1.JRS.9.096031
    wlr=[900 800 700];
    idx=zeros(1,3);
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100*2,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('NIR')
    set(gca,'fontsize',14)
    
    ha(5)=subplot(175);
    % Butz, C., Grosjean, M., Fischer, D., Wunderle, S., Tylmann, W., & Rein, B. (2015). Hyperspectral imaging spectroscopy: a promising method for the biogeochemical analysis of lake sediments. Journal of Applied Remote Sensing, 9, 1–20. https://doi.org/10.1117/1.JRS.9.096031
    wlr=[860 650 555];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100*2,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('CIR')
    set(gca,'fontsize',14)

    ha(6)=subplot(176);
    % Butz, C., Grosjean, M., Fischer, D., Wunderle, S., Tylmann, W., & Rein, B. (2015). Hyperspectral imaging spectroscopy: a promising method for the biogeochemical analysis of lake sediments. Journal of Applied Remote Sensing, 9, 1–20. https://doi.org/10.1117/1.JRS.9.096031
    wlr=[640 545 460];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100*2,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('pRGB')
    set(gca,'fontsize',14)

    ha(7)=subplot(177);
    % Thiele, S. T., Lorenz, S., Kirsch, M., Cecilia Contreras Acosta, I., Tusa, L., Herrmann, E., Möckel, R., & Gloaguen, R. (2021). Multi-scale, multi-sensor data integration for automated 3-D geological mapping. Ore Geology Reviews, 104252. https://doi.org/10.1016/j.oregeorev.2021.104252
    wlr=[680 550 500];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100*2,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('pRGB')
    set(gca,'fontsize',14)

    linkaxes(ha,'xy')
    clear ha
end

if sum(double(strcmp('Images (pRGB, hypRGB, HC1, HC2, L*a*b*, HSV)',index_selec)))
    %% Images
    figure
    ha(1)=subplot(171);
    % Speta, M., Gingras, M. K., & Rivard, B. (2016). Shortwave Infrared Hyperspectral Imaging: A Novel Method For Enhancing the Visibility of Sedimentary And Biogenic Features In Oil-Saturated Core. Journal of Sedimentary Research, 86(7), 830–842. https://doi.org/10.2110/jsr.2016.54
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90)) % 2162 nm = kaolinite, 2199 nm = clay, and 2349 nm = clay
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    title('RGB')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(172);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(rgb2hsv(RGB),-90))
    xlabel('Width (cm)')
    title('HSV')
    set(gca,'fontsize',14)
    
    ha(3)=subplot(173);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(rgb2lab(RGB)/100,-90))
    xlabel('Width (cm)')
    title('LAB')
    set(gca,'fontsize',14)
    
    ha(4)=subplot(174);
    wlr=[1600 1700 2250];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('hypRGB')
    set(gca,'fontsize',14)
    
    ha(5)=subplot(175);
    % Scafutto, R. D. P. M., de Souza Filho, C. R., & Rivard, B. (2016). Characterization of mineral substrates impregnated with crude oils using proximal infrared hyperspectral imaging. Remote Sensing of Environment, 179, 116–130. https://doi.org/10.1016/j.rse.2016.03.033
    wlr=[1722 1760 2311];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('HC1')
    set(gca,'fontsize',14)
    
    ha(6)=subplot(176);
    % Scafutto, R. D. P. M., de Souza Filho, C. R., & Rivard, B. (2016). Characterization of mineral substrates impregnated with crude oils using proximal infrared hyperspectral imaging. Remote Sensing of Environment, 179, 116–130. https://doi.org/10.1016/j.rse.2016.03.033
    wlr=[1722 2311 2349];
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('HC2')
    set(gca,'fontsize',14)

    ha(7)=subplot(177);
    %Thiele, S. T., Lorenz, S., Kirsch, M., Cecilia Contreras Acosta, I., Tusa, L., Herrmann, E., Möckel, R., & Gloaguen, R. (2021). Multi-scale, multi-sensor data integration for automated 3-D geological mapping. Ore Geology Reviews, 104252. https://doi.org/10.1016/j.oregeorev.2021.104252
    wlr=[2200 2250 2350]; % Clay Carbonates Micas
    for i=1:3
        [~,a]=find(abs(wlr(i)-wl)==min(abs(wlr(i)-wl)));
        idx(i)=a;
    end
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(M(:,:,idx)/100,-90))
    colormap(jet)
    xlabel('Width (cm)')
    title('pRGB')
    set(gca,'fontsize',14)

    linkaxes(ha,'xy')
    clear ha
end

%% Aberrant pixels
mask=AbberantPixels(M,RGB,d,0);
maskd=reshape(mask,[],1);
assignin('base','mask',mask);
if sum(double(strcmp('Aberrant pixels',index_selec)))
    figure;
    ha(1)=subplot(121);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(122);
    imagesc(d(1:size(M,1))-d(1),d,imrotate(mask,-90))
    clrmap=[1 0 0;0 1 0];
    colormap(clrmap)
    colorbar('XTick',[0 1],'XTickLabel',{'Aberrant','Relevant'})
    caxis([0 1])
    set(gca,'fontsize',14)
    xlabel('Width (cm)')
    title('Aberrant pixels')
    linkaxes(ha,'xy')
    clear ha
end

if sum(double(strcmp('SNR (local)',index_selec)))
    %% SNR local
    snr=zeros(round(size(RGB,1)/5)-1,round(size(RGB,2)/5)-1);
    for i=1:round(size(RGB,1)/5)-1
        for j=1:round(size(RGB,2)/5)-1
            m=median(reshape(rgb2gray(RGB((i-1)*5+1:(i-1)*5+5,(j-1)*5+1:(j-1)*5+5,:)),[],1));
            s=std(reshape(rgb2gray(RGB((i-1)*5+1:(i-1)*5+5,(j-1)*5+1:(j-1)*5+5,:)),[],1));
            snr(i,j)=m/s;
        end
    end
    SNR=imresize(snr,[size(M,1),size(M,2)]);
    
    figure;
    ha(1)=subplot(121);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(122);
    imagesc(d(1:size(M,1))-d(1),d,imrotate(SNR,-90))
    colormap(jet)
    colorbar
    xlabel('Width (cm)')
    title('SNR (local)')
    set(gca,'fontsize',14)
    linkaxes(ha,'xy')
    clear ha
end

if sum(double(strcmp('Median and standard deviation spectra',index_selec)))
    %% Median and standard deviation spectra
    figure
    yyaxis left
    plot(wl,median(squeeze(M(round(size(M,1)),:,:))),'linewidth',2)
    ylabel('Reflectance (%)')
    yyaxis right
    plot(wl,std(squeeze(M(round(size(M,1)),:,:))),'linewidth',2)
    xlabel('Wavelength (nm)')
    ylabel('Reflectance (%)')
    legend('Median spectrum','Standard deviation spectrum')
    title('Median and standard deviation spectra')
    grid on
    set(gca,'fontsize',14)
end

if sum(double(strcmp('10 most different spectra (Kennard and Stone)',index_selec)))
    %% '10 most different spectra (Kennard and Stone)'
    Md=squeeze(M(round(size(M,1)/2),:,:));
    if mean(wl)<1000
        wlr=[550 975];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
    else
        wlr=[1000 2500];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
    end
    [model,~,~]=kenstone(Md(~isnan(squeeze(mask(round(size(M,1)/2),:))),ind1:ind2),10);
    
    figure
    subplot(211)
    plot(wl,Md(model,:),'linewidth',2)
    ylabel('Reflectance (%)')
    xlabel('Wavelength (nm)')
    title('10 most different spectra (Kennard and Stone)')
    grid on
    xlim([wl(1) wl(end)])
    set(gca,'fontsize',14)
    subplot(212)
    [~,a]=find(mean(Md(model,ind1:ind2),2)>85);
    model(a)=[];
    [~,a]=find(mean(Md(model,ind1:ind2),2)<15);
    model(a)=[];
    plot(wl(ind1:ind2),continuum_removal(wl(ind1:ind2),Md(model,ind1:ind2)),'linewidth',2)
    ylabel('Normalized reflectance')
    xlabel('Wavelength (nm)')
    title('10 most different spectra (Kennard and Stone)')
    grid on
    xlim([wl(1) wl(end)])
    set(gca,'fontsize',12)
    
    clear Md
end

if sum(double(strcmp('Grayscale histogram',index_selec)))
    %% Grayscale histogram   
    figure;
    histogram(reshape(rgb2gray(RGB),[],1),round(length(reshape(rgb2gray(RGB),[],1))/1000),'EdgeColor',[0 0.4470 0.7410])
    xlabel('Greyscale')
    ylabel('Frequency')
    set(gca,'fontsize',14)
end

if sum(double(strcmp('Grayscale histogram as a function of depth',index_selec)))
    %% Grayscale histogram as a function of depth
    figure;
    num=floor(size(RGB,2)/9);
    for i=1:9
        subplot(3,3,i)
        histogram(reshape(rgb2gray(RGB(:,1+(i-1)*num:i*num,:)),[],1),round(length(reshape(rgb2gray(RGB(:,1+(i-1)*num:i*num,:)),[],1))/1000),'EdgeColor',[0 0.4470 0.7410])
        xlim([0 1])
        xlabel('Grayscale')
        ylabel('Frequency')
        grid on
        title([num2str(round(d(1+(i-1)*num),2)),'-',num2str(round(d(i*num),2)),'cm'])
        set(gca,'fontsize',14)
    end
end

if sum(double(strcmp('Spectral histogram',index_selec)))
    %% Spectral histogram
    c=0:.1:100;
    h=zeros(size(M,3),length(c));
    for i=1:size(M,3)
        h(i,:)=hist(squeeze(M(round(size(M,1)/2),:,i)),c);
    end
    h=h';
    hb=h(:);
    
    limi=max([0 median(hb(hb>0))-3*std(hb(hb>0))]);
    limf=min([100 median(hb(hb>0))+3*std(hb(hb>0))]);
    
    [a,~]=find(sum(h,2)>0);
    
    figure;
    imagesc(wl,c(1:max(a)),h(1:max(a),:))
%     surf(wl,c(1:max(a)),h(1:max(a),:))
    caxis([limi limf])
    xlabel('Wavelength (nm)')
    ylabel('Reflectance (%)')
    colormap(jet)
    colorbar
    set(gca,'fontsize',14,'ydir','normal')
end

if sum(double(strcmp('Spectral histogram as a function of depth',index_selec)))
    %% Spectral histogram as a function of depth
    figure;
    c=0:0.1:100;
    num=floor(size(RGB,2)/9);
    ha=[];
    for i=1:9
        ha(1)=subplot(3,3,i);
        h=zeros(size(M,3),length(c));
        for j=1:size(M,3)
            h(j,:)=hist(reshape(M(round(size(M,1)/5*2):round(size(M,1)/5*3),1+(i-1)*num:i*num,j),[],1),c);
        end
        h=h';
        hb=h(:);
        
        limi=max([0 median(hb(hb>0))-4*std(hb(hb>0))]);
        limf=min([100 median(hb(hb>0))+4*std(hb(hb>0))]);
        
        [a,~]=find(sum(h,2)>0);
        
        imagesc(wl,c(1:max(a)),h(1:max(a),:))
        %     surf(wl,c(1:max(a)),h(1:max(a),:))
%         caxis([limi limf])
        xlabel('Wavelength (nm)')
        ylabel('Reflectance (%)')
        title([num2str(round(d(1+(i-1)*num),2)),'-',num2str(round(d(i*num),2)),'cm'])
        colormap(jet)
        colorbar
        set(gca,'fontsize',14,'ydir','normal')
    end
    linkaxes(ha,'xy')
end

if sum(double(strcmp('Raw spectra (2D)',index_selec)))
    %% Raw spectra
    figure;
    ha(1)=subplot(131);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(1,3,2:3);
    imagesc(wl,d,squeeze(median(M,1)))
    title('Raw spectra')
    colormap(flipud(jet))
    h=colorbar;
    h.Label.String = "Reflectance (%)";
    clear h
    xlabel('Wavelength (nm)')
    linkaxes(ha,'y')
    clear ha
    set(gca,'fontsize',14)
end

if sum(double(strcmp('Raw spectra (3D)',index_selec)))
    figure;
    hb(1)=subplot(131);
    imagesc(d(1:size(RGB,1))-d(1),d(ia),imrotate(RGB(:,ia,:),-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    hb(2)=subplot(1,3,2:3);
    surf(wl,d(ia),squeeze(median(M(:,ia,:),1)),'EdgeColor','none')
    title('Raw spectra')
    colormap(flipud(jet))
    h=colorbar;
    h.Label.String = "Reflectance (%)";
    clear h
    ylabel('Depth (cm)')
    xlabel('Wavelength (nm)')
    zlabel('Reflectance (%)')
    linkaxes(hb,'y')
    clear hb
    set(gca,'fontsize',14,'Ydir','reverse')
end

if sum(double(strcmp('Continuum removed (2D)',index_selec)))||...
        sum(double(strcmp('Continuum removed (2D enhanced)',index_selec)))||...
        sum(double(strcmp('Continuum removed (3D)',index_selec)))||...
        sum(double(strcmp('Continuum (2D)',index_selec)))||...
        sum(double(strcmp('Continuum (3D)',index_selec)))
    %% Continuum removed
    [Scr,cr,~]=continuum_removal(wl,squeeze(median(M,1)));
    assignin('base','Scr',Scr);
    assignin('base','cr',cr);
    
    if sum(double(strcmp('Continuum removed (2D)',index_selec)))
        figure;
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(1,3,2:3);
        imagesc(wl,d,Scr),
        title('Continuum removed spectra')
        colormap(flipud(jet))
        h=colorbar;
        h.Label.String = "Normalized reflectance";
        clear h
        caxis([nanmedian(Scr(:))-3*nanstd(Scr(:)) 1])
        xlabel('Wavelength (nm)')
        linkaxes(ha,'y')
        clear ha
        set(gca,'fontsize',14)
    end
    
    if sum(double(strcmp('Continuum removed (2D enhanced)',index_selec)))
        figure;
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(1,3,2:3);
        imagesc(wl,d,Scr),
        title('Continuum removed spectra')
        colormap(flipud(jet))
        h=colorbar;
        h.Label.String = "Normalized reflectance";
        clear h
        caxis([0.8 1])
        xlabel('Wavelength (nm)')
        linkaxes(ha,'y')
        clear ha
        set(gca,'fontsize',14)
    end
    
    if sum(double(strcmp('Continuum removed (3D)',index_selec)))
        figure;
        hb(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d(ia),imrotate(RGB(:,ia,:),-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        hb(2)=subplot(1,3,2:3);
        surf(wl,d(ia),Scr(ia,:),'EdgeColor','none')
        title('Continuum removed spectra')
        colormap(flipud(jet))
        h=colorbar;
        h.Label.String = "Normalized reflectance";
        clear h
        ylabel('Depth (cm)')
        xlabel('Wavelength (nm)')
        zlabel('Reflectance (%)')
        linkaxes(hb,'y')
        clear hb
        set(gca,'fontsize',14,'Ydir','reverse')
    end
    
    %% Continuum
    if sum(double(strcmp('Continuum (2D)',index_selec)))
        figure;
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(1,3,2:3);
        imagesc(wl,d,cr),
        title('Continuum spectra')
        colormap(flipud(jet))
        h=colorbar;
        h.Label.String = "Reflectance (%)";
        clear h
        caxis([nanmedian(cr(:))-3*nanstd(cr(:)) nanmedian(cr(:))+3*nanstd(cr(:))])
        xlabel('Wavelength (nm)')
        linkaxes(ha,'y')
        clear ha
        set(gca,'fontsize',14)
    end
    
    if sum(double(strcmp('Continuum (3D)',index_selec)))
        figure;
        hb(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d(ia),imrotate(RGB(:,ia,:),-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        hb(2)=subplot(1,3,2:3);
        surf(wl,d(ia),cr(ia,:),'EdgeColor','none')
        title('Continuum spectra')
        colormap(flipud(jet))
        h=colorbar;
        h.Label.String = "Reflectance (%)";
        clear h
        ylabel('Depth (cm)')
        xlabel('Wavelength (nm)')
        zlabel('Reflectance (%)')
        linkaxes(hb,'y')
        clear hb
        set(gca,'fontsize',14,'Ydir','reverse')
    end
end

if sum(double(strcmp('FDS (2D)',index_selec)))||...
        sum(double(strcmp('FDS (3D)',index_selec)))
    %% FDS
    Sfds=savgol(squeeze(median(M,1)),7,2,1);
    assignin('base','Sfds',Sfds);
    
    if sum(double(strcmp('FDS (2D)',index_selec)))
        figure;
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(1,3,2:3);
        imagesc(wl,d,Sfds)
        title('First derivative spectra')
        colormap(flipud(jet))
        colorbar
        caxis([nanmedian(Sfds(:))-2*nanstd(Sfds(:)) nanmedian(Sfds(:))+2*nanstd(Sfds(:))])
        xlabel('Wavelength (nm)')
        linkaxes(ha,'y')
        clear ha
        set(gca,'fontsize',14)
    end
    
    if sum(double(strcmp('FDS (3D)',index_selec)))
        figure;
        hb(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d(ia),imrotate(RGB(:,ia,:),-90))
        xlabel('Width (cm)')
        ylabel('Depth (cm)')
        set(gca,'fontsize',14)
        
        hb(2)=subplot(1,3,2:3);
        surf(wl,d(ia),Sfds(ia,:),'EdgeColor','none')
        title('First derivative spectra')
        colormap(flipud(jet))
        colorbar
        ylabel('Depth (cm)')
        xlabel('Wavelength (nm)')
        zlabel('Reflectance (%)')
        linkaxes(hb,'y')
        clear hb
        set(gca,'fontsize',14,'Ydir','reverse')
    end
end

if sum(double(strcmp('Wavelength correlation',index_selec)))
    %% Wavelength correlation
    figure
    imagesc(wl,wl,corr(reshape(M,[],size(M,3))))
    colormap(jet)
    h=colorbar;
    h.Label.String = "Correlation";
    clear h
    ylabel('Wavelength (nm)')
    xlabel('Wavelength (nm)')
    title('Wavelength correlation')
    set(gca,'fontsize',14)
end

if sum(double(strcmp('Q700/500 vs L* (global)',index_selec)))||...
        sum(double(strcmp('Q700/500 vs L* (central)',index_selec)))||...
        sum(double(strcmp('Q700/500 vs L* (map)',index_selec)))
    [RGB,R700_500,L]=Spec2rgb2Q75L(M/100,d,wl);
    assignin('base','RGB',RGB);
    assignin('base','R700_500',R700_500);
    assignin('base','L',L);
    if sum(double(strcmp('Q700/500 vs L* (global)',index_selec)))
        figure
        plot(L.*mask(2:end-1,2:end-1),R700_500.*mask(2:end-1,2:end-1),'b.')
        grid on
        xlabel('L* (%)')
        ylabel('Q700/500')
        title('Q700/500 vs L* (global)')
        set(gca,'fontsize',14)
    end
    if sum(double(strcmp('Q700/500 vs L* (central, top=blue, base=red)',index_selec)))
        figure;
        clr=jet(size(L,2));
        for j=1:size(L,2)
            plot(L(round(size(L,1)/2),j).*mask(round(size(L,1)/2),j+1),R700_500(round(size(L,1)/2),j).*mask(round(size(L,1)/2),j+1),'.','color',clr(j,:))
            hold on
        end
        grid on
        xlabel('L* (%)')
        ylabel('Q700/500')
        title('Q700/500 vs L* (central)')
        set(gca,'fontsize',14)
    end
    if sum(double(strcmp('Q700/500 vs L* (map)',index_selec)))
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(R700_500.*mask(2:end-1,2:end-1),-90))
        h=colorbar;
        h.Label.String = "Q700/500";
        clear h
        caxis([nanmedian(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))-...
            3*nanstd(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))...
            nanmedian(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))+...
            3*nanstd(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))])
        title('700/500 ratio')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(L.*mask(2:end-1,2:end-1),-90))
        caxis([max([0 nanmedian(reshape(L.*mask(2:end-1,2:end-1),[],1))-...
            3*nanstd(reshape(L.*mask(2:end-1,2:end-1),[],1))])...
            min([nanmedian(reshape(L.*mask(2:end-1,2:end-1),[],1))+...
            3*nanstd(reshape(L.*mask(2:end-1,2:end-1),[],1)) 100])])
        colormap(jet)
        h=colorbar;
        h.Label.String = "L*";
        clear h
        title('L* (%)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        clear ha
    end
end

if sum(double(strcmp('Q700/500 vs L* with classification',index_selec)))
    [RGB,R700_500,L]=Spec2rgb2Q75L(M/100,d,wl,1,1);
    assignin('base','RGB',RGB);
    assignin('base','R700_500',R700_500);
    assignin('base','L',L);
end

if sum(double(strcmp('Principal Component Analysis (PCA)',index_selec)))
    if mean(wl)<1000
        %% VNIR PCA D2 pigments
        CompDisp=1:2;
        wlr_vnir=[25 98-5]; % pigments 600-950 => A REVOIR!!!
        
        Sc=savgol(squeeze(M(round(size(M,1)/2),mask(round(size(M,1)/2),:)==1,wlr_vnir(1):wlr_vnir(2))),13,2,2);
        Sd=reshape(M,[],size(M,3));
        S=savgol(Sd(:,wlr_vnir(1):wlr_vnir(2)),13,2,2);
        
        [coeff,~,~,~,explained,~] = pca(Sc(:,8:end-7),'NumComponents',5);
        
        figure;
        subplot(121)
        plot(wl(wlr_vnir(1)+7:wlr_vnir(2)-7),coeff(:,CompDisp),'linewidth',2)
        xlim([wl(wlr_vnir(1)+7), wl(wlr_vnir(2)-7)])
        grid on
        title('VNIR D2')
        legend(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'),...
            strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Wavelength (nm)')
        ylabel('Loading values')
        set(gca,'fontsize',14)
        
        tmp3=S-mean(Sc);
        tmp4=tmp3(:,8:end-7)*coeff(:,CompDisp);
        
        subplot(122)
        biplot(coeff(:,CompDisp),'VarLabels',num2str(wl(wlr_vnir(1)+7:wlr_vnir(2)-7)'))%,'scores',tmp4
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        set(gca,'fontsize',14)
        
        figure
        scatter(tmp4(:,1).*maskd,tmp4(:,2).*maskd,[],'.')
        hold on
        plot(nanmedian(tmp4(:,1).*maskd),nanmedian(tmp4(:,2).*maskd),'w.','markersize',20)
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        title('Scores of the PCA')
        set(gca,'fontsize',14)
        
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,1),size(M,1),size(M,2)).*mask,-90))
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(1)), ' scores'];
        clear h
        caxis([nanmedian(tmp4(:,1).*maskd)-3*nanstd(tmp4(:,1).*maskd) ...
            nanmedian(tmp4(:,1).*maskd)+3*nanstd(tmp4(:,1).*maskd)])
        title(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,2),size(M,1),size(M,2)).*mask,-90))
        caxis([nanmedian(tmp4(:,2).*maskd)-3*nanstd(tmp4(:,2).*maskd) ...
            nanmedian(tmp4(:,2).*maskd)+3*nanstd(tmp4(:,2).*maskd)])
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(2)), ' scores'];
        clear h
        title(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        clear ha
    else
        %% SWIR PCA D1 mineral
        CompDisp=[2 3];
        wlr_swir=[102 125]; %=> A REVOIR!!!
        
        Sc=savgol(squeeze(M(round(size(M,1)/2),mask(round(size(M,1)/2),:)==1,wlr_swir(1):wlr_swir(2))),7,2,1);
        Sd=reshape(M,[],size(M,3));
        S=savgol(Sd(:,wlr_swir(1):wlr_swir(2)),7,2,1);
        
        [coeff,~,~,~,explained,~] = pca(Sc(:,5:end-4),'NumComponents',5);
        
        figure;
        subplot(121)
        plot(wl(wlr_swir(1)+4:wlr_swir(2)-4),coeff(:,CompDisp),'linewidth',2)
        xlim([wl(wlr_swir(1)+4), wl(wlr_swir(2)-4)])
        grid on
        title('SWIR D1')
        legend(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'),...
            strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Wavelength (nm)')
        ylabel('Loading values')
        set(gca,'fontsize',14)
        
        tmp3=S-mean(Sc);
        tmp4=tmp3(:,5:end-4)*coeff(:,CompDisp);
        
        subplot(122)
        biplot(coeff(:,CompDisp),'VarLabels',num2str(wl(wlr_swir(1)+4:wlr_swir(2)-4)'))%,'scores',tmp4
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        set(gca,'fontsize',14)
        
        figure
        scatter(tmp4(:,1).*maskd,tmp4(:,2).*maskd,[],'.')
        hold on
        plot(nanmedian(tmp4(:,1).*maskd),nanmedian(tmp4(:,2).*maskd),'w.','markersize',20)
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        title('Scores of the PCA')
        set(gca,'fontsize',14)
        
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,1).*maskd,size(M,1),size(M,2)),-90))
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(1)), ' scores'];
        clear h
        caxis([nanmedian(tmp4(:,1).*maskd)-3*nanstd(tmp4(:,1).*maskd) ...
            nanmedian(tmp4(:,1).*maskd)+3*nanstd(tmp4(:,1).*maskd)])
        title(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,2).*maskd,size(M,1),size(M,2)),-90))
        caxis([nanmedian(tmp4(:,2).*maskd)-3*nanstd(tmp4(:,2).*maskd) ...
            nanmedian(tmp4(:,2).*maskd)+3*nanstd(tmp4(:,2).*maskd)])
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(2)), ' scores'];
        clear h
        title(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        clear ha
        
        %% SWIR PCA D1 mineral large
        CompDisp=[2 3];
        wlr_swir=[102 144-10]; % => A REVOIR!!!
        
        Sc=savgol(squeeze(M(round(size(M,1)/2),mask(round(size(M,1)/2),:)==1,wlr_swir(1):wlr_swir(2))),7,2,1);
        Sd=reshape(M,[],size(M,3));
        S=savgol(Sd(:,wlr_swir(1):wlr_swir(2)),7,2,1);
        
        [coeff,~,~,~,explained,~] = pca(Sc(:,5:end-4),'NumComponents',5);
        
        figure;
        subplot(121)
        plot(wl(wlr_swir(1)+4:wlr_swir(2)-4),coeff(:,CompDisp),'linewidth',2)
        xlim([wl(wlr_swir(1)+4), wl(wlr_swir(2)-4)])
        grid on
        title('SWIR D1')
        legend(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'),...
            strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Wavelength (nm)')
        ylabel('Loading values')
        set(gca,'fontsize',14)
        
        tmp3=S-mean(Sc);
        tmp4=tmp3(:,5:end-4)*coeff(:,CompDisp);
        
        subplot(122)
        biplot(coeff(:,CompDisp),'VarLabels',num2str(wl(wlr_swir(1)+4:wlr_swir(2)-4)'))%,'scores',tmp4
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        set(gca,'fontsize',14)
        
        figure
        scatter(tmp4(:,1).*maskd,tmp4(:,2).*maskd,[],'.')
        hold on
        plot(nanmedian(tmp4(:,1).*maskd),nanmedian(tmp4(:,2).*maskd),'w.','markersize',20)
        grid on
        xlabel(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        ylabel(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        title('Scores of the PCA')
        set(gca,'fontsize',14)
        
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,1).*maskd,size(M,1),size(M,2)),-90))
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(1)), ' scores'];
        clear h
        caxis([nanmedian(tmp4(:,1).*maskd)-3*nanstd(tmp4(:,1).*maskd) ...
            nanmedian(tmp4(:,1).*maskd)+3*nanstd(tmp4(:,1).*maskd)])
        title(strcat('PC',num2str(CompDisp(1)),':',num2str(round(explained(CompDisp(1)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(tmp4(:,2).*maskd,size(M,1),size(M,2)),-90))
        caxis([nanmedian(tmp4(:,2).*maskd)-3*nanstd(tmp4(:,2).*maskd) ...
            nanmedian(tmp4(:,2).*maskd)+3*nanstd(tmp4(:,2).*maskd)])
        colormap(jet)
        h=colorbar;
        h.Label.String = ['PC',num2str(CompDisp(2)), ' scores'];
        clear h
        title(strcat('PC',num2str(CompDisp(2)),':',num2str(round(explained(CompDisp(2)))),'%'))
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        clear ha
    end
end

if sum(double(strcmp('Minimum Noise Fraction (MNF) + Superpixels',index_selec)))
    %% MNF
    [~, A, noiseFractions] = hyperMnf(squeeze(M(round(size(M,1)/2),mask(round(size(M,1)/2),:)==1,:))',size(M,2),1);
    
    mnf=(reshape(M,[],98)*A)';
    
    figure;
    subplot(121)
    plot(wl,A(:,end:-1:end-2),'linewidth',2)
    xlim([wl(1), wl(end)])
    grid on
    title('MNF components')
    legend(strcat('Noise MNF',num2str(1),':',num2str(round(noiseFractions(end),2))),...
        strcat('Noise MNF',num2str(2),':',num2str(round(noiseFractions(end-1),2))),...
        strcat('Noise MNF',num2str(3),':',num2str(round(noiseFractions(end-2),2))))
    xlabel('Wavelength (nm)')
    ylabel('Loading values')
    set(gca,'fontsize',14)
    
    subplot(122)
    biplot(A(:,end:-1:end-1),'VarLabels',num2str(wl'))%,'scores',tmp4
    grid on
    xlabel(strcat('Noise MNF',num2str(1),':',num2str(round(noiseFractions(end),2))))
    ylabel(strcat('Noise MNF',num2str(2),':',num2str(round(noiseFractions(end-1),2))))
    set(gca,'fontsize',14)
    
    figure
    scatter(mnf(end,:).*maskd,mnf(end-1,:).*maskd,[],'.')
    hold on
    plot(nanmedian(mnf(end,:).*maskd),nanmedian(mnf(end-1,:).*maskd),'w.','markersize',20)
    grid on
    xlabel(strcat('Noise MNF',num2str(1),':',num2str(round(noiseFractions(end),2))))
    ylabel(strcat('Noise MNF',num2str(2),':',num2str(round(noiseFractions(end-1),2))))
    title('Scores of the PCA')
    set(gca,'fontsize',14)
    
    figure
    ha(1)=subplot(131);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    ylabel('Depth (cm)')
    xlabel('Width (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(132);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(mnf(end,:).*maskd,size(M,1),size(M,2)),-90))
    colormap(jet)
    h=colorbar;
    h.Label.String = strcat('Noise MNF',num2str(1),':',num2str(round(noiseFractions(end),2)));
    clear h
    caxis([nanmedian(mnf(end,:).*maskd)-3*nanstd(mnf(end,:).*maskd) ...
        nanmedian(mnf(end,:).*maskd)+3*nanstd(mnf(end,:).*maskd)])
    title(strcat('Noise MNF',num2str(1),':',num2str(round(noiseFractions(end),2))))
    xlabel('Width (cm)')
    set(gca,'fontsize',14)
    
    ha(3)=subplot(133);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(reshape(mnf(end-1,:).*maskd,size(M,1),size(M,2)),-90))
    caxis([nanmedian(mnf(end-1,:).*maskd)-3*nanstd(mnf(end-1,:).*maskd) ...
        nanmedian(mnf(end-1,:).*maskd)+3*nanstd(mnf(end-1,:).*maskd)])
    colormap(jet)
    h=colorbar;
    h.Label.String = strcat('Noise MNF',num2str(2),':',num2str(round(noiseFractions(end-1),2)));
    clear h
    title(strcat('Noise MNF',num2str(2),':',num2str(round(noiseFractions(end-1),2))))
    xlabel('Width (cm)')
    set(gca,'fontsize',14)
    linkaxes(ha,'xy')
    clear ha
    
    mnf=reshape(mnf',size(M,1),size(M,2),size(M,3));
    assignin('base','mnf',mnf);
    assignin('base','A',A);
    assignin('base','noiseFractions',noiseFractions);
    
    mnf=mnf(:,:,end-2:end);
    mnf_min=squeeze(min(min(mnf)));
    mnf_max=squeeze(max(max(mnf)));
    mnf_min_m=ones(size(M,1),size(M,2),3);
    mnf_max_m=ones(size(M,1),size(M,2),3);
    for i=1:3
        mnf_min_m(:,:,i)=squeeze(mnf_min_m(:,:,i))*mnf_min(i);
        mnf_max_m(:,:,i)=squeeze(mnf_max_m(:,:,i))*mnf_max(i);
    end
    MNF=(mnf-mnf_min_m)./(mnf_max_m-mnf_min_m);
    
    [L,~] = superpixels(MNF(:,:,1:3),500,'Compactness',1,'NumIterations',100);
    BW = boundarymask(L);
    
    figure;
    ha(1)=subplot(1,2,1);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(1,2,2);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(imoverlay(RGB,BW,'cyan'),-90))
    xlabel('Width (cm)')
    title('Superpixel')
    set(gca,'fontsize',14)
    linkaxes(ha,'xy')
    clear ha
end

if sum(double(strcmp('Clustering (Kmeans, HAC)',index_selec)))
    %% Clustering
    % S
    evaS = evalclusters(squeeze(median(M,1)),'kmeans','CalinskiHarabasz','KList',1:10);
    % Scr
    [Scr,cr,~]=continuum_removal(wl,squeeze(median(M,1)));
    evaScr = evalclusters(Scr,'kmeans','CalinskiHarabasz','KList',1:10);
    % Cr
    evaCr = evalclusters(cr,'kmeans','CalinskiHarabasz','KList',1:10);
    % Sfds
    Sfds=savgol(squeeze(median(M,1)),7,2,1);
    evaSfds = evalclusters(Sfds,'kmeans','CalinskiHarabasz','KList',1:10);
    
    Km=[evaS.OptimalY evaScr.OptimalY evaCr.OptimalY evaSfds.OptimalY];
    Kmxlab={strcat('S (',num2str(evaS.OptimalK),')'),...
        strcat('Scr (',num2str(evaScr.OptimalK),')'),...
        strcat('Cr (',num2str(evaCr.OptimalK),')'),...
        strcat('Sfds (',num2str(evaSfds.OptimalK),')')};
    assignin('base','Clust_Kmeans',Km);
    
    % S
    evaS = evalclusters(squeeze(median(M,1)),'linkage','CalinskiHarabasz','KList',1:10);
    % Scr
    evaScr = evalclusters(Scr,'linkage','CalinskiHarabasz','KList',1:10);
    % Cr
    evaCr = evalclusters(cr,'linkage','CalinskiHarabasz','KList',1:10);
    % Sfds
    evaSfds = evalclusters(Sfds,'linkage','CalinskiHarabasz','KList',1:10);
    
    HACm=[evaS.OptimalY evaScr.OptimalY evaCr.OptimalY evaSfds.OptimalY];
    HACmxlab={strcat('S (',num2str(evaS.OptimalK),')'),...
        strcat('Scr (',num2str(evaScr.OptimalK),')'),...
        strcat('Cr (',num2str(evaCr.OptimalK),')'),...
        strcat('Sfds (',num2str(evaSfds.OptimalK),')')};
    assignin('base','Clust_HAC',HACm);
    
    figure
    ha(1)=subplot(131);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(132);
    imagesc(1:4,d,Km)
    colormap(ha(2),jet(max(Km(:))))
    h=colorbar;
    h.Label.String = "Number of class";
    clear h
    caxis([0.5 max(Km(:))+0.5])
    title('Kmeans')
    set(gca,'fontsize',14,'xtick',1:4,'xticklabel',Kmxlab)
    
    ha(3)=subplot(133);
    imagesc(1:4,d,HACm)
    colormap(ha(3),jet(max(HACm(:))))
    h=colorbar;
    h.Label.String = "Number of class";
    clear h
    caxis([0.5 max(HACm(:))+0.5])
    title('Agglomerative Hierarchical Tree')
    set(gca,'fontsize',14,'xtick',1:4,'xticklabel',HACmxlab)
    linkaxes(ha,'y')
    clear ha
end

if sum(double(strcmp('Texture analysis',index_selec)))
    %% Texture analysis
    [~,SI] = graycomatrix(rgb2gray(RGB));
    assignin('base','Texture',SI);
    
    figure;
    ha(1)=subplot(121);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(RGB,-90))
    colorbar
    xlabel('Width (cm)')
    ylabel('Depth (cm)')
    set(gca,'fontsize',14)
    
    ha(2)=subplot(122);
    imagesc(d(1:size(RGB,1))-d(1),d,imrotate(SI,-90)),
    colormap(jet),
    h=colorbar;
    h.Label.String = "Co-occurence level";
    clear h
    xlabel('Width (cm)')
    title('Texture analysis with gray-level co-occurrence matrix')
    set(gca,'fontsize',14)
    linkaxes(ha,'xy')
    clear ha
end

end