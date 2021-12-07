function pred=PLS(S,Y, fig, ical, ival)
% Function to create PLSR

% INPUT:
%       S: Spectra dataset
%       Y: Reference dataset

% Set definition
if nargin<4
    [ical,ival]=SetDefinition(size(S,1),round(1/4*size(S,1)));
end
Scal=S(ical,:);
Sval=S(ival,:);
Y_cal=Y(ical,:);
Y_val=Y(ival,:);

% PLSR modeling
Sc_cal=Center(Scal);
Sc_val=Centerval(Sval,Scal);

Y_calb=Y_cal;
Y_valb=Y_val;

Y_cal=Center(Y_cal);

[B,Wstar,T,P,Q,W,R2X,R2Y,RMSE]=plsnipals(Sc_cal,Y_cal,0);
R2X=1-R2X;
R2Y=1-R2Y;

% CALIBRATION
[ac, ~]=size(Y_calb);
Y0 = Sc_cal*Wstar*Q;
Y0 = Uncenterval( Y0, Y_calb );

SEC=sqrt(sum((Y0-Y_calb).^2)/ac);
[~, s2]=size(Y_calb);
z1=zeros(1,s2);
for i=1:s2
    z1(:,i)=sum((Y_calb(:,i)-mean(Y_calb(:,i))).^2);
end
z2=sum((Y0-Y_calb).^2);
SDc=std(Y_calb);
R2c=1-(z2./z1);
RPDc=1./sqrt(1-R2c);

% PREDICTION
[av, ~]=size(Y_valb);
Y1 = Sc_val*Wstar*Q;
Y1 = Uncenterval( Y1, Y_calb );

SEP=sqrt(sum((Y1-Y_valb).^2)/av);
[~, s2]=size(Y_valb);
z1=zeros(1,s2);
if size(Y_valb,1)>1
    for i=1:s2
        z1(:,i)=sum((Y_valb(:,i)-mean(Y_valb(:,i))).^2);
    end
else
    for i=1:s2
        z1(:,i)=sum((Y_valb(:,i)-mean(Y_calb(:,i))).^2);
    end
end
[n, ~]=size(Y1);
z2=sum((Y1-Y_valb).^2);
R2p=1-(z2./z1);
SDp=std(Y_valb);
Bias=sum(Y1-Y_valb)/n;
RMSEP=sqrt(SEP.^2+Bias.^2);
RSD=z2/(n-2);
RPDp=SDp./SEP;
[~, yi]=size(Y_valb);
Intercept=zeros(1,yi);
Slope=zeros(1,yi);
for i=1:yi
    mdl = fitlm(Y_valb(:,yi),Y1(:,yi));
    warning('off','all')
    Coeff=table2array(mdl.Coefficients);
    Intercept(1,i)=round(Coeff(1,1));
    Slope(1,i)=Coeff(2,1);
end

% GRAPHICAL OUTPUT
if fig>0
    figure(fig);
    subplot(2,1,1)
    plot(Y_calb,Y0,'.')
    hold on
    plot([0 max(Y_calb)],[0 max(Y_calb)],'k--')
    title('Calibration')
    grid on
    xlabel('Observed')
    ylabel('Predicted')
    subplot(2,1,2)
    plot(Y_valb,Y1,'.')
    hold on
    plot([0 max(Y_valb)],[0 max(Y_valb)],'k--')
    title('Prediction')
    grid on
    xlabel('Observed')
    ylabel('Predicted')
end

% OUTPUT
pred.coeff.SEC=SEC;
pred.coeff.SEP=SEP;
pred.coeff.RMSEP=RMSEP;
pred.coeff.R2c=R2c;
pred.coeff.R2p=R2p;
pred.coeff.SDc=SDc;
pred.coeff.SDp=SDp;
pred.coeff.RSD=RSD;
pred.coeff.RPDc=RPDc;
pred.coeff.RPDp=RPDp;
pred.coeff.Bias=Bias;
pred.coeff.Intercept=Intercept;
pred.coeff.Slope=Slope;
[~, m]=size(T);
pred.coeff.Recap={'nY' length([Y_calb; Y_valb]); 'Ymin' min([Y_calb; Y_valb]); 'Ymax' max([Y_calb; Y_valb]); 'Ymoy' mean([Y_calb; Y_valb]); 'Ymed' median([Y_calb; Y_valb]); 'Ystd' std([Y_calb; Y_valb]); 'nYcal' length(Y_calb); 'Ycalmin' min(Y_calb); 'Ycalmax' max(Y_calb); 'Ycalmoy' mean(Y_calb); 'Ycalmed' median(Y_calb); 'Ycalstd' std(Y_calb); 'nYval' length(Y_valb); 'Yvalmin' min(Y_valb); 'Yvalmax' max(Y_valb); 'Yvalmoy' mean(Y_valb); 'Yvalmed' median(Y_valb); 'Yvalstd' std(Y_valb); 'VL' m; 'R²c' R2c; 'RPDc' RPDc; 'SEC' SEC; 'SDc' SDc; 'R²p' R2p; 'RPDp' RPDp; 'SEP' SEP; 'RMSEP' RMSEP; 'SDp' SDp; 'RSD' RSD;  'Bias' Bias; 'Intercept' Intercept; 'Slope' Slope};
pred.para.T=T;
pred.para.P=P;
pred.para.W=W;
pred.para.Wstar=Wstar;
pred.para.Q=Q;
pred.para.B=B;
pred.para.R2X=R2X;
pred.para.R2Y=R2Y;
pred.para.RMSE=RMSE;
pred.Pred.Scal=Scal;
pred.Pred.Sval=Sval;
pred.Pred.Ycal=Y_calb;
pred.Pred.Ycalpred=Y0;
pred.Pred.Yval=Y_valb;
pred.Pred.Yvalpred=Y1;
pred.para.nbcomp=m;
end