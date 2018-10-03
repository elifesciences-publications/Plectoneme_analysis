%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single-molecule DNA supercoil data analyzer
%%%
%%% by SHKIM 
%%% Last update 03OCT2018
%%%
%%% Refer to Ganji and Kim et al., 2016 Nano Letters, 16, 4699.
%%% based on a script written by Marijn van Loenhout and by Jacob Kerssemakers 2013  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all;


%% USER parameters------------------------------------------------------------
% Change following parameters before running the codes
initval.Px2um=0.1190391;   % actual image size per pixel. i.e. 1px = 0.119 um
initval.N_Cy5_image=10;   % number of frames taken under red excitation
DNAinfo.DNAlen_bp=20452;  % DNA length in bp
initval.SecondsPerFrame=0.1; % Frame time of the movie (exposure time) in sec.
%% End of USER prameters------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER parameters 

%% Image file settings
DirectoryPath=input('Tiff file path.. ','s');
if isempty(DirectoryPath)
    DirectoryPath=pwd;
end

%% Parameter settings
initval.plecdir=DirectoryPath;
Nframe2read_coil=0;  % set 0 for all frames
Nframe2read_nick=0;  % set 0 for all frames
skip_image=100;  % skip images update during kymograph building 

% DNA box 
initval.boxwidth=7;
initval.boxlen=65;  % must be even number
initval.BGbox_size=5;

% Plot setting
initval.pt_density_imscale=1;   % Unit: kb

initval.interpolation_ratio=100;  % number of oversampleing for interpolation
MG_auto_run=0;  % this is for data re-analysis purpose.

%% prepare matlab
warning off;
release_name=version('-release'); % check MATLAB version


%% --------------------------------------------------------------------------

%% Build_filelist
sdirectory=initval.plecdir;
filelist = dir([initval.plecdir '/*.tiff']);
if isempty(filelist)
    filelist = dir([initval.plecdir '/*.tif']);
end
[Nimages,~]=size(filelist);
disp([num2str(Nimages) ' tif files found']);


%% read first 100 frames for both Cy5 ans SxO
for fr_id=1:initval.N_Cy5_image
    TifLink = Tiff([initval.plecdir '\' filelist(fr_id).name] , 'r');
    Im=double(TifLink.read());
    Im=Im-median(Im(:));
    Im=Im./max(Im(:));
    if fr_id==1
        ImCy5=Im;
    else
        ImCy5=ImCy5+Im;
    end
end
ImCy5=ImCy5/initval.N_Cy5_image;

for fr_id=initval.N_Cy5_image+1:initval.N_Cy5_image*2+1
    TifLink = Tiff([initval.plecdir '\' filelist(fr_id).name] , 'r');
    Im=double(TifLink.read());
    Im=Im-median(Im(:));
    Im=Im./max(Im(:));
    if fr_id==initval.N_Cy5_image+1
        ImCy3=Im;
    else
        ImCy3=ImCy3+Im;
    end
end
ImCy3=ImCy3/initval.N_Cy5_image;



%% make RGB image
Imrgb(:,:,1)=ImCy5+0.2; % 0.08 to make the background bright
Imrgb(:,:,2)=ImCy3+0.2;
Imrgb(:,:,3)=0.2;

% Cy5 to SxO mapping, enable following line if you have bead images for mapping
% Imrgb=m_map_cy5(ImCy3,ImCy5,initval)*3+0.08;    % 0.08 to make the background bright

figure(12)
subplot('Position',[0 0 0.8 1]);
imshow(Imrgb);
imwrite(Imrgb,[initval.plecdir '\merged_rgb.png']);

%% Make an extended cy3 images (original+margin) used for image rotation
[x_size,y_size]=size(ImCy3);
x_size_ext=x_size+(initval.boxwidth+initval.boxlen)*2;
y_size_ext=y_size+(initval.boxwidth+initval.boxlen)*2;
Im_ext = zeros(x_size_ext,y_size_ext);
ext_margin=initval.boxwidth+initval.boxlen;
Im_ext(ext_margin+1:ext_margin+x_size,ext_margin+1:ext_margin+y_size) = ImCy3;

%% load previously analyzed DNA positions
fname_DNApos=[initval.plecdir '\DNAposistions.txt'];
if exist(fname_DNApos,'file')==2
    DNA_positions=load(fname_DNApos);
    [N_DNA_analyzed,~]=size(DNA_positions);
else
    N_DNA_analyzed=0;
    DNA_positions=[];
    % for those analized in the old days...
    filelist_kymo = dir([initval.plecdir '/kymo_*']);
    for fki=1:length(filelist_kymo)
        N_DNA_analyzed=N_DNA_analyzed+1;
        load([initval.plecdir '\' filelist_kymo(fki).name '\allvariables.mat'],'-mat','startpt','endpt','tiltangle','leftend','rightend');
        if exist('leftend','var')
            startpt=leftend(2:-1:1);
            endpt=rightend(2:-1:1);
        end
        tmp_output=[startpt(1) startpt(2) endpt(1) endpt(2) tiltangle];
        DNA_positions=[DNA_positions; tmp_output];
    end
end

% draw the analyzed DNA positions
figure(12)
subplot('Position',[0 0 0.8 1]);    % the RGB image is already drawn here
hold on;
for ndi=1:N_DNA_analyzed
    line([DNA_positions(ndi,1),DNA_positions(ndi,3)],[DNA_positions(ndi,2),DNA_positions(ndi,4)],'Color','y');
    plot(DNA_positions(ndi,1),DNA_positions(ndi,2),'oy');
    text(DNA_positions(ndi,1)+5,DNA_positions(ndi,2),num2str(ndi),'Color','y');
end
hold off;

%% select the DNA for analysis
DNA_id_to_ana=input('choose DNA number to analyze (0 for new) ');
if isempty(DNA_id_to_ana),    DNA_id_to_ana=0;  end

if DNA_id_to_ana==-1
    clear all;
    fclose all;
    return
end

if DNA_id_to_ana==0 % new DNA molecule
    %% Get teter position
    startpt=[nan nan];  % cordinate for Cy5 handle position
    endpt=[nan nan];    % cordinate for non-Cy5 handle position
    tiltangle=0;
    
    disp('Select new DNA');
    disp('Mouse button left: DNA end point with Cy5, right: end point wo Cy5, middle: decision made');
    while 1
        % update the image
        figure(12)
        subplot('Position',[0 0 0.8 1]);
        imshow(Imrgb);
        % indicate the position of the selected DNA
        line([startpt(1),endpt(1)],[startpt(2),endpt(2)],'Color','r');
        hold on;
        plot(startpt(1),startpt(2),'oy');
        % indicate the position of the previously analyzed DNA
        for ndi=1:N_DNA_analyzed
            line([DNA_positions(ndi,1),DNA_positions(ndi,3)],[DNA_positions(ndi,2),DNA_positions(ndi,4)],'Color','y');
            plot(DNA_positions(ndi,1),DNA_positions(ndi,2),'oy');
        end
        hold off;


        %% draw the selected DNA (rotation corrected)
        if sum([isnan(startpt) isnan(endpt)]) == 0  % draw only when a DNA is selected
            % define image area and cordinates for rotation
            left_most=startpt(2)-initval.boxlen-initval.boxwidth+ext_margin;
            right_most=startpt(2)+initval.boxlen+initval.boxwidth+ext_margin;
            up_most=startpt(1)-initval.boxlen-initval.boxwidth+ext_margin;
            bottom_most=startpt(1)+initval.boxlen+initval.boxwidth+ext_margin;

            selectedIm=Im_ext(left_most:right_most,up_most:bottom_most);    % selected area of the DNA + Margin
            selectedImRt = imrotate(selectedIm,-tiltangle,'bilinear','crop');   % image rotation
            selectedMol = selectedImRt(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
                initval.boxlen+1:initval.boxlen+1+initval.boxwidth*2);  % selected area for the DNA (rotated)
            % show the DNA
            figure(12)
            subplot('Position',[0.82 0.1 0.17 0.8]);
            imagesc(selectedMol);axis image;axis off;
            title('Circle is up');
        end

        %% Let user do cliking
        [x,y,flag] = ginput(1);        
        % check if the user cliked outside of the image
        if x<1,        x=1;    elseif x>512,        x=512;    end
        if y<1,        y=1;    elseif y>512,        y=512;    end        
        if flag==2  % user satisfied
            break;
        elseif flag==1  % cy5 handle position update
            startpt = round([x,y]);
        elseif flag==3  % non-cy5 handle position update
            endpt = round([x,y]);
        end

        %% calculate the angle
        tmp=endpt-startpt;
        L3=sqrt(tmp(2)^2+tmp(1)^2);
        tiltangle=acos(tmp(2)/L3);
        tiltangle2=asin(tmp(1)/L3);
        tiltangle=radtodeg(tiltangle*sign(tiltangle2));
        if tmp(1)==0 && tmp(2)<0
            tiltangle=180;
        end        
    end
        
    %% save the DNA position
    tmp_output=[startpt(1) startpt(2) endpt(1) endpt(2) tiltangle];
    DNA_positions=[DNA_positions; tmp_output];
    save([initval.plecdir '\DNAposistions.txt'],'DNA_positions','-ascii')

else   % for the case that the user selected an analyzed DNA again
    % prepare the selected DNA image (rotation corrected)
    startpt=DNA_positions(DNA_id_to_ana,1:2);
    endpt=DNA_positions(DNA_id_to_ana,3:4);
    tiltangle=DNA_positions(DNA_id_to_ana,5);
    
    left_most=startpt(2)-initval.boxlen-initval.boxwidth+ext_margin;
    right_most=startpt(2)+initval.boxlen+initval.boxwidth+ext_margin;
    up_most=startpt(1)-initval.boxlen-initval.boxwidth+ext_margin;
    bottom_most=startpt(1)+initval.boxlen+initval.boxwidth+ext_margin;
    
    selectedIm=Im_ext(left_most:right_most,up_most:bottom_most);    
    selectedImRt = imrotate(selectedIm,-tiltangle,'bilinear','crop');    
    selectedMol = selectedImRt(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
        initval.boxlen+1:initval.boxlen+1+initval.boxwidth*2);
    
    % show the selected DNA
    figure(12)
    subplot('Position',[0.82 0.1 0.17 0.8]);
    imagesc(selectedMol);axis image;axis off;
    title('Circle is up');
end


%% build Coiled kymograph
KymoRaw=[]; % kymograph: without any correction of background, bleaching or whatever
KymoBg=[];  % kymograph: Background corrected
if Nframe2read_coil==0  % zero means all frames
    Nframe2read_coil=length(filelist);
end
% read each frames for SxO and extract intensity profile
fhd13=figure(13);
for fr_id=initval.N_Cy5_image+1:Nframe2read_coil
%     if rem(fr_id,21)~=1
        filenameI0 = [initval.plecdir '\' filelist(fr_id).name];
        TifLink = Tiff(filenameI0 , 'r');
        Im=double(TifLink.read());
        Im_ext = zeros(x_size_ext,y_size_ext);
        Im_ext(ext_margin+1:ext_margin+x_size,ext_margin+1:ext_margin+y_size) = Im;
        Im=Im_ext(left_most:right_most,up_most:bottom_most);
        Im=imrotate(Im,-tiltangle,'bilinear','crop');
        subIm = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
            initval.boxlen+1:initval.boxlen+1+initval.boxwidth*2);
        subImBG1 = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
            initval.boxlen+1-initval.boxwidth:initval.boxlen+1+initval.boxwidth);
        subImBG2 = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
            initval.boxlen+1+initval.boxwidth:initval.boxlen+1+initval.boxwidth*3);
        
        %% -------------------------------------------------------------------
        %Collect summed pixels one vector per image (thus pixel counts of DNA)
        sumIm=sum(subIm,2);
        KymoRaw=[KymoRaw sumIm];
        BGmed=median([subImBG1(:); subImBG2(:)]);
        sumImBg=sum(subIm-BGmed,2);
        KymoBg=[KymoBg sumImBg];
        
        %-----------------plot images-----------------------------------------
        if rem(fr_id,skip_image)==0
            if str2num(release_name(1:4)) >= 2014   % it does not work in earlier version
                set(groot,'CurrentFigure',fhd13);   % this is to update image silently in the background
            else
                figure(fhd13);
            end
            imagesc(subIm');axis('equal');
            title([filelist(fr_id).name ' (' num2str(fr_id) '/' num2str(Nframe2read_coil-initval.N_Cy5_image) ')']);
            axis image;
            drawnow
        end
%     end
end

%% Coiled region selection from the kymograph
% Show kymographs
figure(14);
subplot(1,3,1);
imagesc(KymoBg');
title('Select region of interest');
ylabel('Time (frame)');

% get the frame numbers for coiled region
str_t=1;
[~, tmax]=size(KymoBg);
fin_t=tmax;
disp('mouse button left: first frame to analyze, right: end frame, middle: decision made. '); 
while (1)
    KymoRawCoil=KymoRaw(:,str_t:fin_t);
    KymoBgCoil=KymoBg(:,str_t:fin_t);
    
    % draw the selected region for coiled DNA kymograph
    figure(14);
    subplot(1,3,2);
    imagesc(KymoBgCoil');
    ylabel('Time (frame)');
    title('Coiled selected');
    
    % Let the user clicking for region selection
    [~,clk_t,btntype]=ginput(1);
    if btntype==1   % start frame
        if clk_t<1,	str_t=1;
        else	str_t=clk_t;
        end
    elseif btntype==3   % end frame
        if clk_t>tmax,	fin_t=tmax;
        else	fin_t=clk_t;
        end
    else   % decision made
        break;
    end
    str_t=floor(str_t); %  start frame
    fin_t=floor(fin_t); % end frame
end


%% build nicked kymograph if provided from other folder
NickInOther=menu('Kymograph of Nicked','Select from current','Look up another folder');

if NickInOther==2   % case that the nicked DNA data is in another folder
    %% get the folder name for nicked DNA
    % find the folder name for the nicked DNA (check the time stamp from the folder name)
    dir_seperators=strfind(initval.plecdir,'\');
    parentdir=initval.plecdir(1:dir_seperators(end));
    dir_date=initval.plecdir(end-12:end-7);
    dir_time=str2num(initval.plecdir(end-5:end));
    folderlist=dir(parentdir);
    
    % find the very next folder from time the time stamp
    flg=0;
    for fdi=1:length(folderlist)
        if ~isempty(strfind(folderlist(fdi).name,dir_date))
            tmppath=folderlist(fdi).name;
            tmp_dir_time=str2num(tmppath(end-5:end));
            if tmp_dir_time > dir_time
                flg=1;
                break
            end
        end
    end
    
    if flg    % when the next folder is found...
        initval.nickdir=[parentdir folderlist(fdi).name];
    else    % when it is unclear which folder is for the nick
        initval.nickdir=parentdir;  % Just show the parent directory and let the user select
    end
    
    % ask to user if the found folder name is correct
    initval.nickdir=uigetdir(initval.nickdir);
    
    % check the tiff file info
	filelistnick = dir([initval.nickdir '/*.tiff']);
    [Nimagesnick,~]=size(filelistnick);
    disp([num2str(Nimagesnick) ' tif files found for nick']);
    
    %% Get_tetherpos from user (the same way above)
    filenameI0 = [initval.nickdir '\' filelistnick(3).name];
    Im = double(imread(filenameI0));
    [x_size,y_size]=size(Im);
    x_size_ext=x_size+(initval.boxwidth+initval.boxlen)*2;
    y_size_ext=y_size+(initval.boxwidth+initval.boxlen)*2;
    Im_ext = zeros(x_size_ext,y_size_ext);
    ext_margin=initval.boxwidth+initval.boxlen;
    Im_ext(ext_margin+1:ext_margin+x_size,ext_margin+1:ext_margin+y_size) = Im;
    
   
    while 1
        % draw the DNA line
        figure(1211)
%         subplot('Position',[0 0.51 0.8 0.49]);
        subplot('Position',[0 0 0.8 1]);
        imagesc(Im);
        line([startpt(1),endpt(1)],[startpt(2),endpt(2)],'Color','r');
        hold on;plot(startpt(1),startpt(2),'oy');hold off;
        %%
        if sum([isnan(startpt) isnan(endpt)]) == 0
            left_most=startpt(2)-initval.boxlen-initval.boxwidth+ext_margin;
            right_most=startpt(2)+initval.boxlen+initval.boxwidth+ext_margin;
            up_most=startpt(1)-initval.boxlen-initval.boxwidth+ext_margin;
            bottom_most=startpt(1)+initval.boxlen+initval.boxwidth+ext_margin;
            
            selectedIm=Im_ext(left_most:right_most,up_most:bottom_most);            
            selectedImRt = imrotate(selectedIm,-tiltangle,'bilinear','crop');
            selectedMol = selectedImRt(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
                initval.boxlen+1:initval.boxlen+1+initval.boxwidth*2);
            subplot('Position',[0.82 0.1 0.17 0.8]);
            imagesc(selectedMol);axis image;axis off;
            title('Circle is up');
        end
        %% Let user do cliking
        [x,y,flag] = ginput(1);
        if x<1,        x=1;    elseif x>512,        x=512;    end
        if y<1,        y=1;    elseif y>256,        y=256;    end
        
        if flag==2
            break;
        elseif flag==1
            startpt = round([x,y]);
        elseif flag==3
            endpt = round([x,y]);
        end
        
        %% calculate the angle
        tmp=endpt-startpt;
        
        L3=sqrt(tmp(2)^2+tmp(1)^2);
        tiltangle=acos(tmp(2)/L3);
        tiltangle2=asin(tmp(1)/L3);
        tiltangle=radtodeg(tiltangle*sign(tiltangle2));
        if tmp(1)==0 && tmp(2)<0
            tiltangle=180;
        end
        
    end
    close(1211);
    %% build the kymograph for nick
    KymoRawNick=[];
    KymoBgNick=[];
    if Nframe2read_nick==0;  % zero means all frames
        Nframe2read_nick=length(filelistnick);
    end
    for fr_id=1:Nframe2read_nick        
            filenameI0 = [initval.nickdir '\' filelistnick(fr_id).name];
            TifLink = Tiff(filenameI0 , 'r');
            Im=double(TifLink.read());
            Im_ext = zeros(x_size_ext,y_size_ext);
            Im_ext(ext_margin+1:ext_margin+x_size,ext_margin+1:ext_margin+y_size) = Im;
            Im=Im_ext(left_most:right_most,up_most:bottom_most);
            Im=imrotate(Im,-tiltangle,'bilinear','crop');
            subIm = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
                initval.boxlen+1:initval.boxlen+1+initval.boxwidth*2);
            subImBG1 = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
                initval.boxlen+1-initval.boxwidth:initval.boxlen+1+initval.boxwidth);
            subImBG2 = Im(initval.boxlen+1:initval.boxlen+1+initval.boxlen+initval.boxwidth,...
                initval.boxlen+1+initval.boxwidth:initval.boxlen+1+initval.boxwidth*3);
            
            %% -------------------------------------------------------------------
            %Collect summed pixels one vector per image (thus pixel counts of DNA)
            sumIm=sum(subIm,2);
            KymoRawNick=[KymoRawNick sumIm];
            BGmed=median([subImBG1(:); subImBG2(:)]);
            sumImBg=sum(subIm-BGmed,2);
            KymoBgNick=[KymoBgNick sumImBg];
            
            %-----------------plot images-----------------------------------------
            if rem(fr_id,skip_image)==0
                if str2num(release_name(1:4)) >= 2014
                    set(groot,'CurrentFigure',fhd13);
                else
                    figure(fhd13);
                end
                imagesc(subIm');axis('equal');
                title([filelist(fr_id).name ' (' num2str(fr_id) '/' num2str(Nframe2read_nick) ')']);
                axis image;
                drawnow
            end
    end
    
end


%% nicked region selection from the kymograph
if NickInOther==2
    % Show kymograph nicked
    KymoBg=KymoBgNick;
    KymoRaw=KymoRawNick;
    figure(15);
    subplot(1,3,1);
    imagesc(KymoBg');
    title('Select the region of interest');
    ylabel('Time (frame)');
end

% get the frame numbers for nicked region
str_t=1;
[~, tmax]=size(KymoBg);
fin_t=tmax;
while (1)
    KymoRawNick=KymoRaw(:,str_t:fin_t);
    KymoBgNick=KymoBg(:,str_t:fin_t);
    
    if NickInOther==2
        figure(15);
    else
        figure(14);
    end
    
    subplot(1,3,3);
    imagesc(KymoBgNick');
    title('Nicked selected');
    
    [~,clk_t,btntype]=ginput(1);
    
    if btntype==1
        if clk_t<1
            str_t=1;
        else
            str_t=clk_t;
        end
    elseif btntype==3
        if clk_t>tmax
            fin_t=tmax;
        else
            fin_t=clk_t;
        end
    else
        break;
    end
    str_t=floor(str_t);
    fin_t=floor(fin_t);
end
close(13)
close(14)
if NickInOther==2
    close(15)
end


%% bleach correction
bl_factor=mean(KymoBgNick);
KymoBgBlNick=KymoBgNick./repmat(bl_factor,size(KymoBgNick,1),1);
bl_factor=mean(KymoBgCoil);
KymoBgBlCoil=KymoBgCoil./repmat(bl_factor,size(KymoBgCoil,1),1);


%% Kymograph plotting
figure(16);

subplot(1,2,1);
imagesc(KymoBgBlCoil');
ylabel('Time (frame)');
xlabel('Position (px)');
title('Coiled Bg Bl');

subplot(1,2,2);
imagesc(KymoBgBlNick');
xlabel('Position (px)');
title('Nicked Bg Bl');

%% call script for kymograph analysis
run('MG_ana.m');
