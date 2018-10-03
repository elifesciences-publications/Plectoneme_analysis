
function kym=T1_GetKymoProps(kymograph,datatype)
%This function obtains properties of a kymograph of plectoneme diffusion.

%input: 
    %kymograph. 
    %'datatype'  may be 'mappeddata' or anything else ('raw'). 

%actions
    %Values are obtained per time index (per
    %'profile'or 'Prf'; some of these are further averaged as 'Common' parameters
    %positions associated with plectonemes are found by considering outliers 
    %('peaks') per profile by an iterative outlier procedure, the other
    %positions are associated with background 'plateaus'. These are first
    %estimates used to get a 'common background intensity'. 

    %Finally, ALL counts above this common level are  summed as 'plectonemic counts',
    %(note that this will include small contributions from any noise excursion
    %along the plateau, but will also include plectonems that were initially
    %not found with the outlier procedure)
    %All values below the median background are considered 'nonplectonemic'
    %ccounts. These two are then used to obtain a ratio

%output: 
    % structure 'kym' = 
    % Indices.NonPlec: first-estimate indices of pixels belonging to plateau
    % Indices.Plec: same, (large) plectonemes   
  %arrays of parameters per profile:
    % PrfIndexStarts: left start of plateau positions
    % PrfIndexStops: same, right
    % PrfLevelMinL: average dark signal left of tether
    % PrfLevelMinR: same, right
    % PrfLevelPlateau: median intensity of non-plectonemic pixel positions
    % PrfContentRatios: plectonemic percentage base on common background  
  %scalars for all profiles averaged
    % CommonDarkAvLeft: average dark signal left of tether
    % CommonDarkAvRight: same, right
    % CommonDarkMinLeft: minimum dark signal left of tether
    % CommonDarkMinRight: same, right
    % CommonDarkNoiseLeft: standard deviation of minimum along time axis, left
    % CommonDarkNoiseRight: same, right
    % CommonFluoTreshold: level two sigmas of noise above 'dark' level
    % CommonIndexStart: median start of tether
    % CommonIndexStop: median stop of tether
    % CommonLevelPlateau: median intensity of non-plectonemic pixel positions


%run the program without input to see an example analysis on a artificially
%created kymograph

%JacobKers 2015----------------------------------------------------------
%% DEMO mode
if nargin<1 %test mode
    close all;
    xz=linspace(0,8*pi,80);
    kx=linspace(2*pi,6*pi,100);
    kymograph=rand(100,100)+3;
    for ii=1:100
    fakeprof=sin(xz)+5*exp(-1/(0.5*pi)*(xz-kx(ii)).^2);
    kymograph(ii,10:89)=kymograph(ii,10:89)+2+fakeprof;
    end
    kymograph=kymograph';
    datatype='oridata';
end


% figure(6);imagesc(kymograph);

%% 1) Obtain properties per Profile---------------------------------------
    [rr,cc]=size(kymograph);
    if 0, figure; end
    for ii=1:cc
        prf=kymograph(:,ii);
        prfc=prf;
        
%         figure(61);
%         plot(prf);
        
        %render 'zero-values, from the mapping procedure,not to be used in outlier analysis  
        if strcmp(datatype,'mappeddata')                  
            sel=find(prfc==0); 
            prfc(sel)=NaN;
        end
        %assuming the nonzero section is continuous, we correct output
        %indices later on
               
        %find plateau section by iterative removal of 'outliers'
        [flag,cleandata]=T1_DetermineOutliers(prfc,2,0.8,'all',0);
        insel=find(flag==1);
        outsel=find(flag==0);
          
        %Get curve intensity properties:
        kym.PrfLevelMinL(ii)=min(prf(1:ceil(length(prf)/2))); %left minimum
        kym.PrfLevelMinR(ii)=min(prf(ceil(length(prf)/2):end)); %right minimum
        kym.PrfLevelPlateau(ii)=median(cleandata);             % plateaulevel    
        kym.PrfIndexStarts(ii)=min(insel);  %start of plateau
        kym.PrfIndexStops(ii)=max(insel);   %stop of plateau
        
        st(ii)=std(cleandata);

        kym.Indices(ii).NonPlec=insel; 
        plecsel=find(outsel>=min(insel)&(outsel<=max(insel))&...
                    prf(outsel)>kym.PrfLevelPlateau(ii));  
        kym.Indices(ii).Plec=outsel(plecsel);              
        
        if nargin<1 %test mode     
            plot(prf); hold on;
            plot(kym.Indices(ii).NonPlec,prf(kym.Indices(ii).NonPlec), 'ro'); hold on;
            plot(kym.Indices(ii).Plec,prf(kym.Indices(ii).Plec), 'ko', 'MarkerFaceColor', 'k'); hold on;
            plot(0*prf+kym.PrfLevelPlateau(ii),'r-');  hold off;
            legend('profile', 'plateau estimate', 'plectoneme estimate', 'local median');
            pause(0.001);
            
        end
    end
    
  %% 2) Obtain some properties along frame time axis --------------------
    kym.CommonIndexStart=median(kym.PrfIndexStarts); %'start line'
    kym.CommonIndexStop=median(kym.PrfIndexStops);   %'stop line'
    kym.CommonDarkMinLeft=median(kym.PrfLevelMinL);
    kym.CommonDarkNoiseLeft=std(kym.PrfLevelMinL);
    kym.CommonDarkAvLeft=kym.CommonDarkMinLeft...
                      +2*kym.CommonDarkNoiseLeft;  
                    %note that we find average from minimum+(95% interval/2)
    kym.CommonDarkMinRight=median(kym.PrfLevelMinR);
    kym.CommonDarkNoiseRight=std(kym.PrfLevelMinR);
    kym.CommonDarkAvRight=kym.CommonDarkMinRight...
                      +2*kym.CommonDarkNoiseRight;  
                    %note that we find average from minimum+(95% interval/2)
    kym.CommonDarkAvAll=(kym.CommonDarkAvRight+kym.CommonDarkAvLeft)/2;             
    kym.CommonLevelPlateau=median(kym.PrfLevelPlateau);
 
    %set treshold: every count above is valid, below is dark 
    if strcmp(datatype,'mappeddata') 
        kym.CommonFluoTreshold=0;  %since this date is already tresholded
    else
       kym.CommonFluoTreshold=(kym.CommonDarkAvLeft+2*kym.CommonDarkNoiseLeft...
      +kym.CommonDarkAvRight+2*kym.CommonDarkNoiseRight)/2;
    end
    kym.PrfContentRatios=T2_GetRatios(kym,kymograph);                                                 
    kym=orderfields(kym); 
    
 %% Plotting menu   
    if nargin<1
        close all;
        subplot(2,2,1); pcolor(kymograph); shading flat, colormap hot;
        subplot(2,2,2);
        plot(kym.PrfLevelPlateau); hold on;
        plot(kym.PrfLevelMinL,'k'); hold on;
        plot(kym.PrfLevelMinR,'r'); hold on;
        legend('plateau', 'left minimum', 'right minimum');
        subplot(2,2,3);
        plot(kym.PrfIndexStarts); hold on;
        plot(kym.PrfIndexStops); hold on;
        plot(0*kym.PrfIndexStarts+kym.CommonIndexStart,'r');
        plot(0*kym.PrfIndexStops+kym.CommonIndexStop, 'r');
        legend('starts', 'stops');
        subplot(2,2,4);
        plot(kym.PrfContentRatios);
        legend('percentage in excursions');
    end
    
    
    
    function ContentRatios=T2_GetRatios(kym,kymograph)
% example: kym = 
%         CommonDarkAvLeft: -0.0926
%        CommonDarkAvRight: -0.1139
%        CommonDarkMinLeft: -0.1810
%       CommonDarkMinRight: -0.2011
%      CommonDarkNoiseLeft: 0.0442
%     CommonDarkNoiseRight: 0.0436
%       CommonFluoTreshold: -0.0154
%         CommonIndexStart: 9
%          CommonIndexStop: 60
%       CommonLevelPlateau: 1.2440
%                  Indices: [1x62 struct]
%           PrfIndexStarts: []
%            PrfIndexStops: []
%             PrfLevelMinL: []
%             PrfLevelMinR: []
%          PrfLevelPlateau: []


[rr,cc]=size(kymograph);
axz=(1:rr);
ContentRatios=zeros(cc,1);
for ii=1:cc
    Prf=kymograph(:,ii);
    plecsel=kym.Indices(ii).Plec;
    nonplecsel=kym.Indices(ii).NonPlec;     
    
    %remove 'dark'level
    PrfCor=Prf-kym.CommonFluoTreshold;  
    sel=find(PrfCor<0); PrfCor(sel)=0;    
    
    %'shave off' excursions
    nonpleclev=kym.CommonLevelPlateau-kym.CommonFluoTreshold; 
    PrfCorNonPlec=PrfCor; PrfCorNonPlec(plecsel)=nonpleclev;
    
    %isolate excursions (potential plectonemes)
    PrfCorPlec=PrfCor-PrfCorNonPlec;  
    
    
    switch 2
        case 1
        %Base content on equating outliers to plectonemes
        ContentRatios(ii)=100*sum(PrfCorPlec)/sum(PrfCor);
        case 2
        %base content on using all 'shaved' counts above median level
        countsel=find(Prf>nonpleclev);
        ExcessCounts=sum(PrfCor(countsel)-nonpleclev);
        ContentRatios(ii)=100*sum(ExcessCounts)/sum(PrfCor);
    end
    
    if 0
    subplot(1,1,1);    
    plot(Prf); hold on;
    plot(axz(nonplecsel),Prf(nonplecsel),'rx'); hold on;
    plot(axz(plecsel),Prf(plecsel),'bo','Markersize',5, 'MarkerFaceColor', 'b'); hold on; 
    plot(0*Prf+kym.PrfLevelPlateau(ii),'r-');  hold on;    
    
    plot(axz(countsel),Prf(countsel),'ko', 'Markersize',10); hold off; 
    title('Levels');
    legend('profile', 'plateau points estimate', 'plectoneme estimate', 'common median plateau',  'above common plateau');
    % subplot(1,2,2); 
    [~]=ginput(1);
    end
end
    
 %%    
function [flag,cleandata]=T1_DetermineOutliers(data,tolerance,sigchange,how,sho);
%this function is meant to find a representative value for a standard
%deviation in a heavily skewed distribution (typically, flat data with
% %peaks). It calculates the standard deviation and average the data;
% Based on these, outliers are determined and excluded for a new calculation
% of average and SD; this is repeated until sigma does not change anymore too much
% . This is repeated until the new sigma does not change much
% %anymore
%output: positions of outliers

%Jacob Kers 2013 and before---------------------------------------------
binz=50;


if nargin<5  %For testing/demo purposes
    close all
    data=JK00_DEMODATA_Peaks;
    tolerance=2;
    sigchange=0.7;
    how='positive';
    sho=1;
    plot(data,'o-');
    binz=20;
end

sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
    sigma_old=sigma;
    selc=find(flag==1);
    data(flag==1); 
    ls=length(selc);
    av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
    sigma=nanstd(data(selc));
    ratio=sigma/sigma_old;
    switch how
        case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
    end
    %plot menu------------------  
    if sho==1
        cleandata=data(selc); 
        hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
        sthst=hist(cleandata,hx);
        bar(hx,sthst);
        title('Histogram');
        dum=ginput(1);
        pause(0.5);     
    end
    %---------------------------- 
    selc=find(flag==1);
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);





