function [RGBim,R700_500, L]=Spec2rgb2Q75L(M,d,wl,fig,clust)
% Function that convert hyperspectral image into RGB, and to L*a*b* to compute a equivalent to the Q7/4 diagram.
% INPUT:
%           M: Hyperspectral datacube (n*m*p)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
%           classif: To make clustering (clust=1)
% OUTPUT:
%           RGBim: RGB image
%           R700_450: ratio of 700 nm and 475 nm
%           L: Luminance
% Debret, M., Sebag, D., Desmet, M., Balsam, W., Copard, Y., Mourier, B., Susperrigui, A.-S., Arnaud, F., Bentaleb, I., Chapron, E., Lallier-Vergès, E., Winiarski, T., 2011. Spectrocolorimetric interpretation of sedimentary dynamics: The new “Q7/4 diagram.�? Earth-Science Rev. 109, 1–19. https://doi.org/10.1016/j.earscirev.2011.07.002

if size(M,3)>1
    S=reshape(M,[],size(M,3));
end
if median(S(:))>100
    S=S/10000;
    M=M/10000;
end
if nargin<4
    fig=0;
end
if nargin<5
    clust=0;
end
RGB=M(:,:,[50 28 13]);

%% Aberrant pixels
mask=AbberantPixels(M,RGB,d,0);
maskd=reshape(mask,[],1);
assignin('base','mask',mask);

% 1931 CIE color matching functions
[x, y, z, ~, ~] = Get_xyz();
% close all

% Interpolation
nmi=x(1,1):0.01:x(end,1);
xi=interp1(x(:,1)',x(:,2)',nmi);
yi=interp1(y(:,1)',y(:,2)',nmi);
zi=interp1(z(:,1)',z(:,2)',nmi);

wlr=round(wl,2);

% Matching between hyperspectral wavelength and 1931 CIE color
xs=zeros(1,length(wlr));
ys=zeros(1,length(wlr));
zs=zeros(1,length(wlr));
for i=1:length(wlr)
    [~,b]=find(abs(wlr(i)-nmi)==min(abs(wlr(i)-nmi)));
    if abs(wlr(i)-nmi(b))<0.1
        xs(i)=xi(b);
        ys(i)=yi(b);
        zs(i)=zi(b);
    else
        xs(i)=0;
        ys(i)=0;
        zs(i)=0;
    end
end

% Wavelength conversion to CIE XYZ
X=S*xs';
Y=S*ys';
Z=S*zs';

% Normalization
XYZ=[X/sum(xs) Y/sum(ys) Z/sum(zs)];

% CIE XYZ conversion to L*a*b* and RGB
LAB=xyz2lab(XYZ);
RGB=xyz2rgb(XYZ);

% Search for wavelengths for the ratio
wl500=find(abs(wlr-500)==min(abs(wlr-500)));
wl700=find(abs(wlr-700)==min(abs(wlr-700)));

% Creation of the L* and the ratio maps
L=reshape(LAB(:,1),size(M,1),size(M,2));
R700_500=reshape(S(:,wl700)./S(:,wl500),size(M,1),size(M,2));
RGBim=reshape(RGB,size(M,1),size(M,2),3);

% Median filtering preprocessing
L=medfilt2(L,[20,3]);
L=L(2:end-1,2:end-1);
R700_500=medfilt2(R700_500,[20,3]);
R700_500=R700_500(2:end-1,2:end-1);
Ld=reshape(L,[],1);
R700_500d=reshape(R700_500,[],1);

RGBim(1,:)=0;
RGBim(end,:)=0;
RGBim(:,1)=0;
RGBim(:,end)=0;

if fig>0 % Figure display
    if clust==1
        f1=figure;
        plot(L.*mask(2:end-1,2:end-1),R700_500.*mask(2:end-1,2:end-1),'b.')
        grid on
        xlabel('L* (%)')
        ylabel('Q700/500')
        title('Q700/500 vs L*')
        set(gca,'fontsize',14)
        
        Cknn=ginput();
        close(f1)
        
        % KNN
        [idxknn,~]=knnsearch(Cknn,[L(:) R700_500(:)],'Distance','mahalanobis');
        idxknn=reshape(idxknn(:,1),size(M,1)-2,size(M,2)-2);
        D=pdist2(Cknn,[L(:) R700_500(:)],'mahalanobis')';
        D=abs((D-min(D(:)))/(max(D(:))-min(D(:)))-1);
        %         D=abs(D./repmat(sum(D,2),1,size(D,2))-1);
        %         D=(D./repmat(sum(D,2),1,size(D,2)));
        KNN_mix=reshape(D,size(M,1)-2,size(M,2)-2,size(Cknn,1));
        
        % Kmeans
        %         evaS = evalclusters([L(round(size(L,1)/2),:);R700_475(round(size(L,1)/2),:)]','kmeans','CalinskiHarabasz','KList',[1:10]);
        [~,Ckmed,~,~,~,~]= kmedoids([L(round(size(L,1)/2),:);R700_500(round(size(L,1)/2),:)]',size(Cknn,1),'Distance','mahalanobis');
        idxkmed = knnsearch(Ckmed,[L(:) R700_500(:)],'Distance','mahalanobis');
        idxkmed=reshape(idxkmed(:,1),size(M,1)-2,size(M,2)-2);
        D=pdist2(Ckmed,[L(:) R700_500(:)],'mahalanobis')';
        D=abs((D-min(D(:)))/(max(D(:))-min(D(:)))-1);
        %         D=abs(D./repmat(sum(D,2),1,size(D,2))-1);
        %         D=(D-min(D(:)))/(max(D(:))-2*min(D(:)));
        Kmed_mix=reshape(D,size(M,1)-2,size(M,2)-2,size(Ckmed,1));
        
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(R700_500.*mask(2:end-1,2:end-1),-90))
        colorbar
        caxis([nanmedian(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))-...
            3*nanstd(reshape(R700_500.*mask(2:end-1,2:end-1),[],1)) ...
            nanmedian(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))+...
            3*nanstd(reshape(R700_500.*mask(2:end-1,2:end-1),[],1))])
        title('700/500 ratio')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(L.*mask(2:end-1,2:end-1),-90))
        caxis([max([0 nanmedian(reshape(L.*mask(2:end-1,2:end-1),[],1))-...
            3*nanstd(reshape(L.*mask(2:end-1,2:end-1),[],1))])...
            min([nanmedian(reshape(L.*mask(2:end-1,2:end-1),[],1))+...
            3*nanstd(reshape(L.*mask(2:end-1,2:end-1),[],1)) 100])])
        colormap(jet)
        colorbar
        title('L* (%)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        clear ha
        
        figure;
        subplot(1,4,1:2)
        c=lines(size(Cknn,1));
        for i=1:size(Cknn,1)
            plot(Ld(idxknn(:)==i),R700_500d(idxknn(:)==i),'.','color',c(i,:))
            hold on
        end
        for i=1:size(c,1)
            plot(Cknn(i,1),Cknn(i,2),'.','markersize',50,'markeredgecolor',c(i,:))
            hold on
        end
        grid on
        xlabel('L* (%)')
        ylabel('700/500 ratio')
        title('Q75 diagram')
        set(gca,'fontsize',14)
        
        ha2(1)=subplot(143);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        title('RGB image')
        set(gca,'fontsize',14)
        
        ha2(4)=subplot(144);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(idxknn,-90))
        colormap(ha2(4),lines(size(Cknn,1)))
        colorbar('southoutside')
        title('KNN classification map')
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14,'YAxisLocation','right')
        linkaxes(ha2,'xy')
        
        subsize1 = get(ha2(1), 'Position');
        subsize2 = get(ha2(4), 'Position');
        set(ha2(1),'Units','normalized', 'Position', [subsize1(1) subsize2(2) subsize1(3) subsize2(4)]);
        clear ha2
        
        numfig=size(KNN_mix,3);
        
        figure
        ha(1)=subplot(1,numfig+1,1);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        title('RGB image')
        set(gca,'fontsize',14)
        
        for i=1:numfig
            ha(i+1)=subplot(1,numfig+1,i+1);
            imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(squeeze(KNN_mix(:,:,i)),-90))
            colormap(jet)
            colorbar('southoutside')
            title(strcat('Class ',num2str(i),' probability (KNN)'))
            xlabel('Width (cm)')
            set(gca,'fontsize',14)
        end
        linkaxes(ha,'xy')
        
        subsize1 = get(ha(1), 'Position');
        subsize2 = get(ha(i+1), 'Position');
        set(ha(1),'Units','normalized', 'Position', [subsize1(1) subsize2(2) subsize1(3) subsize2(4)]);
        clear ha
        
        for i=1:size(KNN_mix,3)
            [a,b]=find(min(min(squeeze(KNN_mix(:,:,i))))==squeeze(KNN_mix(:,:,i)));
            Snear(i,:)=squeeze(M(a(1),b(1),:));
            Sneard1(i,:)=savgol(squeeze(M(a(1),b(1),:))',7,2,1);
            [a,~]=find(idxkmed(:)==i);
            Smed(i,:)=median(S(a,:));
            Smedd1(i,:)=savgol(median(S(a,:)),7,2,1);
        end
        
        figure
        subplot(221)
        plot(wl(10:end-10),Snear(:,10:end-10),'linewidth',4)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('Reflectance (%)')
        title('Spectra closest to the poles (KNN)')
        set(gca,'fontsize',14)
        
        subplot(222)
        plot(wl(10:end-10),Smed(:,10:end-10),'linewidth',4)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('Reflectance (%)')
        title('Median spectra for each medoids (K Medoids)')
        set(gca,'fontsize',14)
        
        subplot(223)
        plot(wl(10:end-10),Sneard1(:,10:end-10),'linewidth',4)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('First derivative spectra')
        title('Spectra closest to the poles (KNN)')
        set(gca,'fontsize',14)
        
        subplot(224)
        plot(wl(10:end-10),Smedd1(:,10:end-10),'linewidth',4)
        grid on
        xlabel('Wavelength (nm)')
        ylabel('First derivative spectra')
        title('Median spectra for each medoids (K Medoids)')
        set(gca,'fontsize',14)
        
        figure;
        subplot(1,4,1:2)
        for i=1:size(Ckmed,1)
            plot(Ld(idxkmed(:)==i),R700_500d(idxkmed(:)==i),'.')
            hold on
        end
        plot(Ckmed(:,1),Ckmed(:,2),'k.','markersize',50)
        grid on
        xlabel('L* (%)')
        ylabel('700/500 ratio')
        title('Q75 diagram')
        set(gca,'fontsize',14)
        
        ha(1)=subplot(143);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        title('RGB image')
        set(gca,'fontsize',14)
        
        ha(4)=subplot(144);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(idxkmed,-90))
        colormap(ha(4),lines(size(Ckmed,1)))
        colorbar('southoutside')
        title('Kmedioids classification map')
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14,'YAxisLocation','right')
        linkaxes(ha,'xy')
        
        subsize1 = get(ha(1), 'Position');
        subsize2 = get(ha(4), 'Position');
        set(ha(1),'Units','normalized', 'Position', [subsize1(1) subsize2(2) subsize1(3) subsize2(4)]);
        clear ha
        
        figure
        ha(1)=subplot(1,numfig+1,1);
        imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        title('RGB image')
        set(gca,'fontsize',14)
        
        for i=1:numfig
            ha(i+1)=subplot(1,numfig+1,i+1);
            imagesc(d(1:size(RGBim,1))-d(1),d,imrotate(squeeze(Kmed_mix(:,:,i)),-90))
            colormap(jet)
            colorbar('southoutside')
            title(strcat('Class ',num2str(i),' probability (K Medoids)'))
            xlabel('Width (cm)')
            set(gca,'fontsize',14)
        end
        linkaxes(ha,'xy')
        
        subsize1 = get(ha(1), 'Position');
        subsize2 = get(ha(i+1), 'Position');
        set(ha(1), 'Units','normalized','Position', [subsize1(1) subsize2(2) subsize1(3) subsize2(4)]);      
        clear ha
    else
        figure;
        plot(LAB(:,1),S(:,wl700)./S(:,wl500),'.')
        grid on
        xlabel('L* (%)')
        ylabel('Q700/500')
        title('Q700/500 vs L* (global)')
        set(gca,'fontsize',14)
        
        figure;
        plot(L(round(size(L,1)/2),:),R700_500(round(size(L,1)/2),:),'.')
        grid on
        xlabel('L* (%)')
        ylabel('Q700/500')
        title('Q700/500 vs L* (central)')
        set(gca,'fontsize',14)
        
        figure
        ha(1)=subplot(131);
        imagesc(d(1:size(RGBim,1)),d,imrotate(RGBim,-90))
        ylabel('Depth (cm)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(2)=subplot(132);
        imagesc(d(1:size(RGBim,1)),d,imrotate(R700_500,-90))
        colorbar
        caxis([nanmedian(R700_500(:))-3*nanstd(R700_500(:)) ...
            nanmedian(R700_500(:))+3*nanstd(R700_500(:))])
        title('700/500 ratio')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        
        ha(3)=subplot(133);
        imagesc(d(1:size(RGBim,1)),d,imrotate(L,-90))
        caxis([max([0 nanmedian((L(:)))-3*nanstd((L(:)))])...
            min([nanmedian((L(:)))+3*nanstd((L(:))) 100])])
        colormap(jet)
        colorbar
        title('L* (%)')
        xlabel('Width (cm)')
        set(gca,'fontsize',14)
        linkaxes(ha,'xy')
        
        subsize1 = get(ha(1), 'Position');
        subsize2 = get(ha(3), 'Position');
        set(ha(1),'Units','normalized', 'Position', [subsize1(1) subsize2(2) subsize1(3) subsize2(4)]);
        clear ha
    end
end
end