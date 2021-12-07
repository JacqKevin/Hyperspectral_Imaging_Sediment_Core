function [Rv_s, p_s, Ss, Yc] = CorrelIndiceMap(C,dC,Y,dY,figi,pas,pas2)
% Correlation between an indice and chemical analysis
% INPUT :
%           C : Index map(s)
%           dC: Depth vector
%           Y: Reference data to compare with the index
%           dY: Depth of the reference sampling
%           optional:
%               figi : Diplay or not figures (1/0)
%               pas: Depth limits of the reference sampling (Max, Min,
%                   Center)
%               pas2: Width of the reference sampling
% OUTPUT :
%           Rv_s : Correlation between the index and the reference
%           p_s : Associated p-value
%           Ss : Standard deviation of the index
%           Yc : Index sampled to the reference resolution

% Zone ref
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
else if strcmp(dataref,'Min')
        for i=1:length(dY)
            dyb(1,i)=dY(i)+1/3*datadref(i);
            dyb(2,i)=dY(i)+2/3*datadref(i);
        end
    else if strcmp(dataref,'Center')
            for i=1:length(dY)
                dyb(1,i)=dY(i)-1/6*datadref(i);
                dyb(2,i)=dY(i)+1/6*datadref(i);
            end
        end
    end
end

% Compare dm and dy
iter=1;
diy=zeros(2,length(dY));
for i=1:length(dY)
    for j=1:2
        [~, a]=find(abs(dC-dyb(j,i))==min(abs(dC-dyb(j,i))));
        if abs(dC(a(1))-dyb(j,i))<0.2
            dis(j,iter)=a(1);
            Ye(iter,:)=Y(i,:);
        else
            diy(j)=1;
        end
    end
    if abs(dC(a(1))-dyb(j,i))<0.5
        iter=iter+1;
    end
end
if dis(1,1)==1
    dis=dis(:,2:end);
    Ye=Ye(2:end,:);
end
dis(dis==0)=1;
dY=dY(diy(1,:)==0);
Y=Y(diy(1,:)==0,:);

% Remove 1/3 of the sides (edge effects)
Mi=C(round(2/5*size(C,1):3/5*size(C,1)),:);
Mia=C(round(2/5*size(C,1):3/5*size(C,1)),:);

for i=1:size(Ye,1)
    Yi(i,1)=nanmedian(nanmedian(Mi(:,dis(1,i):dis(2,i))));
end

[Rv_s,p_s]=corr(Ye(isnan(Yi)==0),Yi(isnan(Yi)==0));
Ss=median(nanstd(C));

if figi>0
    figure;
    subplot(211)
    if median(C(:))+3*std(C(:))>median(Y)&median(C(:))-3*std(C(:))<median(Y)
    else
        yyaxis left
    end
    plot(dC,nanmedian(C),'linewidth',2)
    ylabel('Index values')
    if median(C(:))+3*std(C(:))>median(Y)&median(C(:))-3*std(C(:))<median(Y)
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
    ylabel('Index values')
    title(strcat('Correlation :',num2str(Rv_s),'(p-value:',num2str((p_s)),'), Standard deviation :',num2str(Ss)))
end

Yc.Y=Ye;
Yc.Yp=Yi;

end