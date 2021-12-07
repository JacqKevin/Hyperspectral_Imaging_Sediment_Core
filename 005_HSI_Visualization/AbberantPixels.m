function mask=AbberantPixels(M,RGB,d,fig,thresh)
% Function to detect abberant pixels based on a threshold in standard
% deviation.
% INPUT:
%           M: Hyperspectral datacube (n*m*p)
%           RGB: Associated RGB image (n*m*3)
%           d: Associated depth (1*m)
%           fig: Display figure (fig>0)
%           thresh: Threshold:
%                  <100: number of standard deviation
%                  >100: reflectance threshold
% OUTPUT:
%           mask: Map of abberant pixels (n*m)

IM=reshape(median(reshape(M,[],size(M,3)),2),size(M,1),size(M,2));
if nargin>5
    if thresh>100
        lim1=thresh;
%         lim2=9000;
    else if thresh<100
            lim1=median(IM(:))-thresh*std(IM(:));
%             lim2=median(IM(:))+lim*std(IM(:));
        end
    end
else
    lim1=median(IM(:))-2*std(IM(:));
    lim2=median(IM(:))+2*std(IM(:));
end
IMd=reshape(IM,[],1);
maskd=ones(size(M,1)*size(M,2),1);
maskd(IMd<lim1|IMd>lim2)=NaN;
maskd(IMd<lim1)=NaN;
mask=reshape(maskd,size(M,1),size(M,2));

if nargin>3
    if fig>0
        figure;
        ha(1)=subplot(311);
        imagesc(d,d(1:size(M,1)),RGB)
        set(gca,'fontsize',14)
        colorbar
        
        ha(2)=subplot(312);
        imagesc(d,d(1:size(M,1)),IM)
        caxis([median(IM(:))-3*std(IM(:)) median(IM(:))+3*std(IM(:))])
        colormap(jet)
        set(gca,'fontsize',14)
        colorbar
        
        ha(3)=subplot(313);
        imagesc(d,d(1:size(M,1)),mask)
        colormap(jet)
        colorbar
        caxis([0 1])
        set(gca,'fontsize',14)
        xlabel('Depth (cm)')
        linkaxes(ha,'xy')
    end
end

end