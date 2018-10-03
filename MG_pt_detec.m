%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script read DNA density kymograph and identify the position and size
% of plectonemes.
% 
% Originally written by Marijn van Loenhout;
% edited/adapted by Jacob Kers 2013
% edited/adapted by SHKIM DEC2013
% Version 8. edited/adapted by SHKIM JUL2016
%
% MG_gen and MG_ana must be run before
% This program uses the result (workspace) from MG_ana
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMPORTANT PARAMETERS
pt.Npx2sum_sizedetermin=2;  % actual size: 1+2*Npx2sum_sizedetermin
pt.size_factor=1.1;    % This is to compensate un-added side tails of gaussian distribution
pt_threshold=0.25;   % actual threshold will be N_base_per_px_Nick*pt.min_pt_size
pt.max_diff_len=3;  % maximum lerteral jump allowed
pt.min_pt_len = 2;

base_resolution=500;    % unit bp
initval.base_resolution=base_resolution;    % (UNIT:bp) binsize for histogram analysis
base_vector=base_resolution/2:base_resolution:DNAinfo.DNAlen_bp;
base_vector=(base_vector/1000)';   % to display in kb
base_vector_2x=(base_resolution/4:base_resolution/2:DNAinfo.DNAlen_bp)'/1000;
base_vector_half=(base_resolution:base_resolution*2:DNAinfo.DNAlen_bp)'/1000;

%% check kymograph and parameters
kymo_pt=Kymo_DNAdensity;
[Nfr_org,Npx_org]=size(kymo_pt);
DNA_start=floor(DNAinfo.start_Coil_avr);
DNA_end=ceil(DNAinfo.start_Coil_avr+DNAinfo.len_Nick_px);

%% draw kymographs
figure(51);clf
figure(51);
subplot(1,3,1);
imagesc(kymo_pt);caxis([0 initval.pt_density_imscale]);
xlabel('Position (pixel)');
ylabel('Time (frame)');
title('DNA density');
drawnow;

pt.first_frame_of_kymo_analyzed=1;
pt.last_frame_of_kymo_analyzed=Nfr_org;
pt.N_frame_of_kymo_analyzed=Nfr_org;

%% smoothed kymograph
% disp('smoothing...(may take some time)');
kymo_pt_sm=Smooth2DimArray(kymo_pt);
figure(51);
subplot(1,3,2);
imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);
title('smoothed');xlabel('position (pixel)');drawnow;

%% get background DNA density and fluctuation

% simple median calculation
kymo_inside=kymo_pt(:,ceil(DNAinfo.start_Coil_avr):floor(DNAinfo.len_Nick_px));
rough_med=median(kymo_inside(:));
BG_DNAdensity_std=std(kymo_inside(:));
backgroundall=kymo_inside((kymo_inside<(rough_med+BG_DNAdensity_std*2)) & kymo_inside>(rough_med-BG_DNAdensity_std*2));
BG_DNAdensity_coil=median(backgroundall);

% Jakob's median calculation
kym=T1_GetKymoProps(Kymo_DNAdensity','mappeddata');
BG_DNAdensity_coil=kym.CommonLevelPlateau;

disp(['estimated background DNA density: ' num2str(BG_DNAdensity_coil)]);

pt.min_pt_size_kb=BG_DNAdensity_coil*pt_threshold;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start plectoneme detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% peakfinding------------------------------------

% disp('start peak finding');
[Nfr,Npx]=size(kymo_pt_sm);
figure(51)
subplot(1,3,3);
imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);
title('peak position');
xlabel('pixel');
drawnow;

kymo_peaks=zeros(Nfr,Npx);
peaks=cell(Nfr,1);

for fri=1:Nfr
    tmp=MG_Get1Dpeaks(kymo_pt_sm(fri,:),BG_DNAdensity_coil,pt_threshold,pt.Npx2sum_sizedetermin);
    
    % remove peaks that is too close to the edge
    tmp = tmp( tmp > DNAinfo.start_Coil_avr);
    tmp = tmp ( tmp < DNAinfo.end_Coil_avr);
    
    peaks(fri)={tmp};
    kymo_peaks(fri,tmp)=1;
end

% draw the peaks found
figure(51); subplot(1,3,3); hold on;
for fri=1:Nfr
    if ~isempty(peaks{fri})
        plot(peaks{fri},fri, 'k.');
    end
end
hold off;drawnow;

%------------------------------------------------------



%% Connect peaks

% disp('connecting peaks');

% put the peaks in frame 1 (first line of kymo) into the new plectoneme array
pts={};
[~,Npts]=size(peaks{1});
pts_start_fr=(1:Npts)*0;
pts_end_fr=(1:Npts)*0;
peaks_ptid=cell(Nfr,1);
for peakid=1:Npts,    % put the peaks found in frame 1 into the new pts
    pts(peakid)={peaks{1}(peakid)};
    peaks_ptid{1}(peakid)=peakid;
    pts_start_fr(peakid)=1;
    pts_end_fr(peakid)=1;
end
pts=pts';   % this is for debuging purpose. Do not delete.
% Run through the frames
for fri=2:Nfr
    cur_peaks=peaks{fri};
    [~,N_cur_peaks]=size(cur_peaks);
    
    for peakid=1:N_cur_peaks
        %1: highlight current peak
%         hold on;
%         plot(cur_peaks(peakid),fri, 'w.');hold off;
        %2: find nearest neighbor in the previous frame

        if isempty(peaks{fri-1})    % check if there is no pt detected in the previous frame
            ptjump=pt.max_diff_len+1;
        else
            distances=abs(cur_peaks(peakid)-peaks{fri-1});
            [ptjump,nearest]=min(distances);            
        end
        
        if ptjump > pt.max_diff_len
            %3: if no nearest neighbor, creast new pt
            Npts=Npts+1;
            pts(Npts)={cur_peaks(peakid)};
            pts_start_fr(Npts)=fri;
            pts_end_fr(Npts)=fri;
            peaks_ptid{fri}(peakid)=Npts;
            pts_peakid(Npts)=peakid;
        else
            %4: if nearest neighbor, connect to the previous pt
            % find which pt it belongs to
            tmp_ptid=peaks_ptid{fri-1}(nearest);
            
            % check if it is connected to another pt
            if pts_end_fr(tmp_ptid)==fri
                % determine which one is more close
                pt1=pts{tmp_ptid}(end);
                pt2=cur_peaks(peakid);
                prev=pts{tmp_ptid}(end-1);
                if abs(prev-pt1) > abs(prev-pt2)
                    % move pt1 to new pts
                    Npts=Npts+1;
                    pts(Npts)={pts{tmp_ptid}(end)};
                    pts_start_fr(Npts)=fri;
                    pts_end_fr(Npts)=fri;
                    pts_peakid(Npts)=pts_peakid(tmp_ptid);
                    peaks_ptid{fri}(pts_peakid(tmp_ptid))=Npts;
                    % connect pt2 to the exist pts
                    pts{tmp_ptid}(end)=cur_peaks(peakid);
                    pts_peakid(tmp_ptid)=peakid;
                    peaks_ptid{fri}(peakid)=tmp_ptid;
                else
                    % make pt2 as new
                    Npts=Npts+1;
                    pts(Npts)={cur_peaks(peakid)};
                    pts_start_fr(Npts)=fri;
                    pts_end_fr(Npts)=fri;
                    pts_peakid(Npts)=peakid;
                    peaks_ptid{fri}(peakid)=Npts;
                end
            else
                % put the peak to exist pts
                pts{tmp_ptid}=[pts{tmp_ptid} cur_peaks(peakid)];
                peaks_ptid{fri}(peakid)=tmp_ptid;
                pts_end_fr(tmp_ptid)=fri;
                pts_peakid(tmp_ptid)=peakid;
            end
        end
    end
end


%% Remove short diffusions
pt.N_pt=0;
clear('tmp_pts','tmppts_start_fr','tmppts_end_fr');

for ttpti=1:Npts
    pt_len=length(pts{ttpti});
    if pt_len >= pt.min_pt_len
        % put into the tmp_pts
        pt.N_pt=pt.N_pt+1;
        tmp_pts(pt.N_pt)=pts(ttpti);
        tmppts_start_fr(pt.N_pt)=pts_start_fr(ttpti);
        tmppts_end_fr(pt.N_pt)=pts_end_fr(ttpti);
    end
end

if pt.N_pt==0
    disp('no plectoneme found');
    return
end

pt.pos_in_px=tmp_pts';
pt.start_fr=tmppts_start_fr;
pts_end_fr=tmppts_end_fr;


%% determine the size and position of the plectoneme in base position
% [tmax,xmax]=size(kymo_pt_sm);
kymo_pts=ones(Nfr,Npx);

pt.pos_in_bp={};
pt.size_in_bp={};
for ttpti=1:pt.N_pt
    c_pt=pt.pos_in_px{ttpti};
    pt_len=length(c_pt);
    for fri=1:pt_len
        tmpx=c_pt(fri);
        tmpt=pt.start_fr(ttpti)+fri-1;
        pt.pos_in_bp{ttpti,1}(fri)=sum(kymo_pt(tmpt,1:tmpx))-kymo_pt(tmpt,tmpx)/2;
        pt.size_in_bp{ttpti,1}(fri)=sum(kymo_pt(tmpt,tmpx-pt.Npx2sum_sizedetermin:tmpx+pt.Npx2sum_sizedetermin)-BG_DNAdensity_coil);
        pt.size_in_bp{ttpti,1}(fri)=pt.size_in_bp{ttpti,1}(fri)*pt.size_factor;
        kymo_pts(tmpt,tmpx)=pt.size_in_bp{ttpti,1}(fri);
    end
end

% analsysis on all the peaks found
peaks_in_bp=peaks;  % This is only for memory allocation
for fri=1:Nfr
    N_peaks_raw=length(peaks{fri});
    for peak_id=1:N_peaks_raw
        peaks_in_bp{fri}(peak_id)=sum(kymo_pt_sm(fri,1:peaks{fri}(peak_id)))-kymo_pt_sm(fri,peaks{fri}(peak_id))/2;
    end
end





%% draw the diffusion
figure(51)
subplot(1,3,3);
% imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);
% title('peak position');xlabel('pixel');

hold on;
multicolormode=0;
if multicolormode
    pt_color=rand(pt.N_pt,3);
    for ttpti=1:pt.N_pt
        pt_len=length(pt.pos_in_px{ttpti});
        timevector=tmppts_start_fr(ttpti)-1+(1:pt_len);
        %     pt_color(ttpti,:)=[rand() rand() rand()];
        p=plot(pt.pos_in_px{ttpti},timevector);
        set(p,'Color',pt_color(ttpti,:),'LineWidth',2);
        %     input('sfjkdljs');
    end
else
    for ttpti=1:pt.N_pt
        pt_len=length(pt.pos_in_px{ttpti});
        timevector=tmppts_start_fr(ttpti)-1+(1:pt_len);
        p=plot(pt.pos_in_px{ttpti},timevector,'r-');
        set(p,'LineWidth',1);
    end
end
hold off;
drawnow;
close(51);



%% --------------------------------------------------------------------------
% save plectonme info and logs

cur_time=clock;
cur_time(6)=floor(cur_time(6));
output_dir=['plec3_' num2str(cur_time(1)) num2str(cur_time(2)) num2str(cur_time(3)) ...
    '_' num2str(cur_time(4)) 'h' num2str(cur_time(5)) 'm' num2str(cur_time(6),2) 's'];

if ~exist('MG_pt_rerun','var')
    mkdir(initval.plecdir,output_dir);
    outfile_name=[initval.plecdir '\' output_dir '\plectoneme.mat'];
else
    outfile_name=[data_path '\' output_dir '\plectoneme.mat'];
    mkdir(data_path,output_dir);
    kymo_output_dir=[data_path '\' kymo_path];
end
if ~exist('kymo_output_dir','var')
    kymo_output_dir='';
end
save(outfile_name,'pt','kymo_output_dir','kymo_pt_sm');



run('MG_pt_play.m');



