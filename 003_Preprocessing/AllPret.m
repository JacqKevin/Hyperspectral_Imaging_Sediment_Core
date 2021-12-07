function Spret = AllPret(S,wl,pret)
% Function to compute several preprocessing algorithms.

% INPUT:    X: Spectra matrix (n*m)
%           wl: Wavelengths vector (1*m)
%           pret(optionnel): To choose a specific preprocessing (0:tous, 1:detrend,
%           2:SNV, 3:SNVD, 4:MSC, 5:D1, 6:D2, 7:SNV+D1, 8:SNVD+D1,
%           9:MSC+D1, 10:SNV+D2, 11:SNVD+D2, 12:MSC+D2, 13:CR)

if nargin<3
    pret=0;
end

if iscell(S)
    S=cell2mat(S);
end

Spret.Raw=S;

if ischar(pret)||iscell(pret)
    if  strcmp(pret,'Detrend')
        pret=1;
    else if strcmp(pret,'SNV')
            pret=2;
        else if strcmp(pret,'SNVD')
                pret=3;
            else if strcmp(pret,'MSC')
                    pret=4;
                else if strcmp(pret,'D1')
                        pret=5;
                    else if strcmp(pret,'D2')
                            pret=6;
                        else if strcmp(pret,'SNV+D1')||strcmp(pret,'SNVD1')
                                pret=7;
                            else if strcmp(pret,'SNVD+D1')||strcmp(pret,'SNVDD1')
                                    pret=8;
                                else if strcmp(pret,'MSC+D1')||strcmp(pret,'MSCD1')
                                        pret=9;
                                    else if strcmp(pret,'SNV+D2')||strcmp(pret,'SNVD2')
                                            pret=10;
                                        else if strcmp(pret,'SNVD+D2')||strcmp(pret,'SNVDD2')
                                                pret=11;
                                            else if strcmp(pret,'MSC+D2')||strcmp(pret,'MSCD2')
                                                    pret=12;
                                                else if strcmp(pret,'CR')
                                                        pret=13;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if pret==0||pret==1
    % Detrend
    if pret>0
        Spret=detrend(S);
    else
        Spret.Detrend=detrend(S);
    end
end

if pret==0||pret==2
    % SNV
    if pret>0
        Spret=snv(S);
    else
        Spret.SNV=snv(S);
    end
end

if pret==0||pret==3
    % SNVD
    if pret>0
        Spret=snv(detrend(S));
    else
        Spret.SNVD=snv(detrend(S));
    end
end

if pret==0||pret==4
    % MSC
    if pret>0
        if size(S,1)>1
            Spret=msc(S);
        else
            Spret=S;
        end
    else
        if size(S,1)>1
            Spret.MSC=msc(S);
        else
            Spret.MSC=S;
        end
    end
end

if pret==0||pret==5
    % D1
    if pret>0
        Spret=savgol(S,7,2,1);
    else
        Spret.D1=savgol(S,7,2,1);
    end
end

if pret==0||pret==6
    % D2
    if pret>0
        Spret=savgol(S,9,2,2);
    else
        Spret.D2=savgol(S,9,2,2);
    end
end

if pret==0||pret==7
    % SNV+D1
    if pret>0
        Spret=savgol(snv(S),7,2,1);
    else
        Spret.SNVD1=savgol(snv(S),7,2,1);
    end
end

if pret==0||pret==8
    % SNVD+D1
    if pret>0
        Spret=savgol(snv(detrend(S)),7,2,1);
    else
        Spret.SNVDD1=savgol(snv(detrend(S)),7,2,1);
    end
end

if pret==0||pret==9
    % MSC+D1
    if pret>0
        if size(S,1)>1
            Spret=savgol(msc(S),7,2,1);
        else
            Spret=savgol((S),7,2,1);
        end
    else
        if size(S,1)>1
            Spret.MSCD1=savgol(msc(S),7,2,1);
        else
            Spret.MSCD1=savgol((S),7,2,1);
        end
    end
end

if pret==0||pret==10
    % SNV+D2
    if pret>0
        Spret=savgol(snv(S),9,2,2);
    else
        Spret.SNVD2=savgol(snv(S),9,2,2);
    end
end

if pret==0||pret==11
    % SNVD+D2
    if pret>0
        Spret=savgol(snv(detrend(S)),9,2,2);
    else
        Spret.SNVDD2=savgol(snv(detrend(S)),9,2,2);
    end
end

if pret==0||pret==12
    % MSC+D2
    if pret>0
        if size(S,1)>1
            Spret=savgol(msc(S),9,2,2);
        else
            Spret=savgol((S),9,2,2);
        end
    else
        if size(S,1)>1
            Spret.MSCD2=savgol(msc(S),9,2,2);
        else
            Spret.MSCD2=savgol((S),9,2,2);
        end
    end
end

if pret==0||pret==13
    % Continuum Removal
    [S_cr,cr]=continuum_removal(wl',S');
    if pret>0
        Spret=S_cr;
    else
        Spret.CR=S_cr;
        Spret.CRc=cr;
    end
end

end