function [Rv_s, p_s, Ss, Yc] = CorrelCartoY(CartoY,dCartoY,Y,dY,figi,pas,pas2)
% Correlation between an the predicted map and chemical analysis
% INPUT :
%           CartoY : Prediction map(s)
%           dCartoY: Depth vector
%           Y: Reference data to compare with the prediction map
%           dY: Depth of the reference sampling
%           optional:
%               figi : Diplay or not figures (1/0)
%               pas: Width of the reference sampling
%               pas2: Depth limits of the reference sampling (Max, Min,
%                   Center)
% OUTPUT :
%           Rv_s : Correlation between the prediction and the reference
%           p_s : Associated p-value
%           Ss : Standard deviation of the prediction
%           Yc : Index sampled to the reference resolution

% Subsampling
if pas==0
    dataref = questdlg('What is the sampling position?','References','Min','Center','Max','Max');
    
    datadref = inputdlg('What is the sampling width in centimeters?');
    datadref=str2double(datadref);
else
    if length(pas)==1
        pas=repmat(pas,1,length(dY));
    end
    datadref=pas;
    dataref=pas2;
end

dyb=zeros(2,length(dY));
if strcmp(dataref,'Max')
    for i=1:length(dY)
        dyb(1,i)=dY(i)-2/3*datadref(i);
        dyb(2,i)=dY(i)-1/3*datadref(i);
    end
else
    if strcmp(dataref,'Min')
        for i=1:length(dY)
            dyb(1,i)=dY(i)+1/3*datadref(i);
            dyb(2,i)=dY(i)+2/3*datadref(i);
        end
    else
        if strcmp(dataref,'Center')
            for i=1:length(dY)
                dyb(1,i)=dY(i)-1/6*datadref(i);
                dyb(2,i)=dY(i)+1/6*datadref(i);
            end
        end
    end
end

% Compare and et dy
iter=1;
diy=zeros(2,length(dY));
for i=1:length(dY)
    for j=1:2
        [~, a]=find(abs(dCartoY-dyb(j,i))==min(abs(dCartoY-dyb(j,i))));
        if abs(dCartoY(a(1))-dyb(j,i))<0.2
            dis(j,iter)=a(1);
            Ye(iter,:)=Y(i,:);
        else
            diy(j)=1;
        end
    end
    if abs(dCartoY(a(1))-dyb(j,i))<0.5
        iter=iter+1;
    end
end
if dis(1,1)==1
    dis=dis(:,2:end);
    Ye=Ye(2:end,:);
end
if dis(1,1)==0
    dis(1,1)=1;
end
dY=dY(diy(1,:)==0);
Y=Y(diy(1,:)==0,:);

% Remove 1/3 of the sides (edge effects)
Mi=CartoY(round(2/5*size(CartoY,1):3/5*size(CartoY,1)),:);

for i=1:size(Ye,1)
    Yi(i,1)=nanmedian(nanmedian(Mi(:,dis(1,i):dis(2,i))));
end

[Rv_s,p_s]=corr(Ye(~isnan(Yi)),Yi(~isnan(Yi)));
Ss=nanmedian(nanstd(CartoY));

if figi>0
    figure;
    subplot(211)
    if median(CartoY(:))+3*std(CartoY(:))>median(Y)&&median(CartoY(:))-3*std(CartoY(:))<median(Y)
    else
        yyaxis left
    end
    plot(dCartoY,nanmedian(CartoY),'linewidth',2)
    ylabel('Predicted')
    if median(CartoY(:))+3*std(CartoY(:))>median(Y)&&median(CartoY(:))-3*std(CartoY(:))<median(Y)
        hold on
    else
        yyaxis right
    end
    plot(dY,Y,'.','markersize',20)
    ylabel('Observed')
    grid on
    subplot(212)
    plot(Ye,Yi,'.','markersize',20)
    grid on
    xlabel('Observed')
    ylabel('Predicted')
    title(strcat('Correlation :',num2str(Rv_s),'(p-value:',num2str((p_s)),'), Standard deviation :',num2str(Ss)))
end

Yc.Y=Ye;
Yc.Yp=Yi;

end