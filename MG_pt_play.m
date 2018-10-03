%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script analyzes kymograph generated from kymo_gen_sh
% Originally written by Marijn van Loenhout;
% edited/adapted by Jacob Kers 2013
% edited/adapted by SHKIM DEC2013
% Version 8 by SHKIM on JUL2016
%
% kymo_gen_sh, kymo_ana_sh and kymo_pt_detec must be run before.
% This program uses the result (workspace) from MG_pt_detec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define local parameters
% base_per_pixel=DNAinfo.N_base_per_px_Nick;
time_unit=initval.SecondsPerFrame; % sec
N_pt=pt.N_pt;


%% Draw kymograph and plectoneme position
figure(71)
subplot(1,2,1);
imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);
ylabel('Time (frame)');
title('DNA density (unit: kb)','FontSize',12);
% colorbar('location','NorthOutside');
xlabel('Position (pixel)');
subplot(1,2,2);
imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);
hold on;
for ttpti=1:pt.N_pt
    pt_len=length(pt.pos_in_px{ttpti});
    timevector=pt.start_fr(ttpti)-1+(1:pt_len);
    p=plot(pt.pos_in_px{ttpti},timevector,'k-','LineWidth',2);
end
hold off;
% colorbar('location','NorthOutside');
title('plectoneme position','FontSize',12);
xlabel('Position (pixel)');



%% Basic plectoneme parameters determined

% calc mean number of pt per unit time
pt.mean_N_pt=0;
for ttpti=1:pt.N_pt
    pt.mean_N_pt=pt.mean_N_pt+length(pt.pos_in_px{ttpti});
end
pt.mean_N_pt=pt.mean_N_pt/pt.N_frame_of_kymo_analyzed;



%% build binary plectoneme occupation kymograph
[Kymo_pt_occupation_half,~]=build_Kymo_pt_occupation(pt,DNAinfo,base_vector_half);
[Kymo_pt_occupation,~]=build_Kymo_pt_occupation(pt,DNAinfo,base_vector);
[Kymo_pt_occupation_2x,~]=build_Kymo_pt_occupation(pt,DNAinfo,base_vector_2x);



%% plotting the plectonem position in pixel space
figure(72);
% plot on top of the kymograph
subplot(1,2,1);
imagesc(kymo_pt_sm);caxis([0 initval.pt_density_imscale]);

pt_color=rand(pt.N_pt,3);
hold on;
for ttpti=1:pt.N_pt
    pt_len=length(pt.pos_in_px{ttpti});
    timevector=pt.start_fr(ttpti)-1+(1:pt_len);
    p=plot(pt.pos_in_px{ttpti},timevector,'k-','LineWidth',2);
    set(p,'Color',pt_color(ttpti,:),'LineWidth',2);
end
hold off;
xlabel('Position (pixel)');
ylabel('Time (frame)');
title('DNA density');

%% plot pt_occupation kymograph and pt positions in basepair space
colormap gray;
subplot(1,2,2,'replace');
% imagesc(base_vector,1:pt.N_frame_of_kymo_analyzed,~Kymo_pt_occupation');
hold on;
for ttpti=1:pt.N_pt
    pt_len=length(pt.pos_in_bp{ttpti});
    timevector=pt.start_fr(ttpti)-1+(1:pt_len);
    p=plot(pt.pos_in_bp{ttpti},timevector,'LineWidth',0.1);
%     p=plot(pt.pos_in_bp{ttpti}*base_per_pixel/1000,timevector,'LineWidth',0.1);
    set(p,'Color',pt_color(ttpti,:),'LineWidth',2);
end
hold off;
xlabel('position in DNA (kbp)');
ylim([1 pt.N_frame_of_kymo_analyzed]);
axis ij;
xlim([0 DNAinfo.DNAlen_bp/1000]);
ylabel('Time (frame)');
box on
title('Plectonemic DNA');

%% tmppp
% 
% tmp_kymo_show=1:150;
% 
% figure(72);
% % plot on top of the kymograph
% subplot(1,2,1);
% imagesc(kymo_pt_sm(tmp_kymo_show,:));caxis([0 initval.pt_density_imscale]);
% pt_color=rand(pt.N_pt,3);
% hold on;
% for ttpti=1:pt.N_pt
%     pt_len=length(pt.pos_in_px{ttpti});
%     timevector=pt.start_fr(ttpti)-1+(1:pt_len);
%     p=plot(pt.pos_in_px{ttpti},timevector,'k-','LineWidth',2);
% %     set(p,'Color',pt_color(ttpti,:),'LineWidth',2);
%     set(p,'Color',[1 0 0],'LineWidth',1);
% end
% hold off;
% xlabel('Position (pixel)');
% ylabel('Time (frame)');
% title('DNA density');
% 
% colormap gray;
% subplot(1,2,2,'replace');
% imagesc(base_vector,1:length(tmp_kymo_show),Kymo_pt_occupation(:,tmp_kymo_show)');
% hold on;
% for ttpti=1:pt.N_pt
%     pt_len=length(pt.pos_in_bp{ttpti});
%     timevector=pt.start_fr(ttpti)-1+(1:pt_len);
%     p=plot(pt.pos_in_bp{ttpti},timevector,'LineWidth',0.1);
% %     p=plot(pt.pos_in_bp{ttpti}*base_per_pixel/1000,timevector,'LineWidth',0.1);
%     set(p,'Color',pt_color(ttpti,:),'LineWidth',1);
%     set(p,'Color',[1 0 0],'LineWidth',1);
% end
% hold off;
% xlabel('position in DNA (kbp)');
% % ylim([1 pt.N_frame_of_kymo_analyzed]);
% axis ij;
% xlim([0 DNAinfo.DNAlen_bp/1000]);
% ylabel('Time (frame)');
% box on
% title('Plectonemic DNA');

%%
% 
% plot(base_vector,mean(Kymo_pt_occupation(:,tmp_kymo_show),2));
% xlabel('position in DNA (kbp)');
% % ylim([1 pt.N_frame_of_kymo_analyzed]);
% xlim([0 DNAinfo.DNAlen_bp/1000]);
% ylabel('Plectoneme occupancy');
% box on
% % title('Plectonemic DNA');
% 
% plot(base_vector,mean(Kymo_pt_occupation(:,tmp_kymo_show),2));
% xlabel('position in DNA (kbp)');
% % ylim([1 pt.N_frame_of_kymo_analyzed]);
% xlim([0 DNAinfo.DNAlen_bp/1000]);
% ylabel('Plectoneme occupancy');
% box on
% title('Plectonemic DNA');


%% position dependent plectoneme lifetime
pt=update_pt_life_time(pt,initval,base_vector);


%% Diffusion Coeeficient
% [pt,msd]=calc_DiffCoef(pt,initval,min_trace_len,min_pt_size_ana)
[pt, msd, dsp]=calc_DiffCoef(pt,initval,10,0.3);

if isfield(pt,'D_mean_px')
%     disp(['Diff. coef. (mean of indv.): ' num2str(pt.D_mean_px,2) '+' num2str(pt.D_se_px,2) ' um^2/s (' num2str(msd.n_valid_MSD) ' pts)']);
%     disp(['Diff. coef. (mean of indv.): ' num2str(pt.D_mean_bp,2) '+' num2str(pt.D_se_bp,2) ' kbp^2/s (' num2str(msd.n_valid_MSD) ' pts)']);
%     disp(['Diff. coef. of averaged MSD : ' num2str(pt.D_of_meanMSD_freeDiff_px,2) ' um^2/s (' num2str(msd.n_valid_MSD) ' pts)']);
%     disp(['Diff. coef. of averaged MSD : ' num2str(pt.D_of_meanMSD_freeDiff_bp,2) ' kbp^2/s (' num2str(msd.n_valid_MSD) ' pts)']);
else
    disp('No msd calculated');
end

% plot displacement and MSD
% in px space
% figure(731);
% subplot(2,2,1);
% plot(dsp.times_easy,dsp.px_easy);
% xlabel('time (s)','fontsize',12);
% ylabel('Displacement (\mum)','fontsize',12);
% subplot(2,2,2);
% plot(msd.times_easy,msd.easy_px);
% xlabel('time (s)','fontsize',12);
% ylabel('MSD (\mum^2/s)','fontsize',12);
% title([int2str(N_pt) ' plectonemes'],'fontsize',12);
% hold off;
% in bp space
figure(732);
subplot(2,2,1);
plot(dsp.times_easy,dsp.bp_easy);
xlabel('time (s)','fontsize',12);
ylabel('Displacement (kbp)','fontsize',12);
subplot(2,2,2);
plot(msd.times_easy,msd.easy_bp);
xlabel('time (s)','fontsize',12);
ylabel('MSD (kbp^2/s)','fontsize',12);
title([int2str(N_pt) ' plectonemes'],'fontsize',12);
hold off;

% plot MSDs that were used to calc D
if isfield(msd,'n_valid_MSD')   % sometimes, if there's no pt traces that are longer enough to calculated MSD, no diffusion coefficient will be calculated
    figure(732)
    subplot(2,2,3);
    plot(msd.analyzed_time,msd.analyzed_bp);
    xlabel('time (s)','fontsize',12);
    ylabel('MSD (kbp^2/s)','fontsize',12);
    title([num2str(msd.n_valid_MSD) ' plctonemes (>' num2str(msd.min_trace_len_ana*time_unit) 's, >' num2str(msd.min_pt_size_ana) 'kb)'],'fontsize',12);
    ylim([0 max(msd.analyzed_bp(:))*1.2]);
    text(msd.tmax*time_unit*0.4,max(msd.analyzed_bp(:)),['D_i_d_v=' num2str(pt.D_mean_bp,2) ' kbp^2/s']);
    
    subplot(2,2,4);
    errorbar((1:length(msd.mean_msd_bp))*time_unit,msd.mean_msd_bp,msd.analyzed_se_bp,'k.-');hold on;
    plot(msd.msd_fit_linear,'r');
    plot(msd.msd_fit_sub_diffusion,'g');hold off;
    ylabel('mean MSD (kbp^2/s)','fontsize',12);
    ylim([0 max(msd.mean_msd_bp)*1.3])
    xlim([0 (msd.min_trace_len_ana+1)*time_unit*1.2])
    xlabel('Time (s)','fontsize',12);
    text(msd.min_trace_len_ana*time_unit*0.5,max(msd.mean_msd_bp)*0.2,['D_a_v_r=' num2str(pt.D_of_meanMSD_freeDiff_bp,2) ' kbp^2/s']);
    text(msd.min_trace_len_ana*time_unit*0.1,max(msd.mean_msd_bp)*1.15,['D_s_u_b=' num2str(pt.D_of_meanMSD_SubDiff_bp,2) ' kbp^2/s']);
    text(msd.min_trace_len_ana*time_unit*0.1,max(msd.mean_msd_bp)*0.95,['n_s_u_b=' num2str(pt.n_of_meanMSD_SubDiff_bp,2)]);
    legend off;    
end

%% calc plectoneme density and occupancy
[pt_prof]=build_ptden(pt,base_vector,Kymo_pt_occupation);
[pt_prof_2x]=build_ptden(pt,base_vector_2x,Kymo_pt_occupation_2x);


%% Size (in base) vs. Position (in base) and pt density
figure(75);
subplot(2,1,1,'replace');
plot(0);    %   There may be more sophistcated way to do this, but this is to draw a solid border of the graph
hold on;
pt_color=rand(N_pt,3);
for pti=1:N_pt
%     p=plot(pt.pos_in_bp{pti},pt.size_in_bp{pti});
    p=plot(pt.pos_in_bp{pti},pt.size_in_bp{pti},'.','MarkerSize',10);
%     set(p,'Color',pt_color(pti,:),'LineWidth',0.5);
    xlim([1 5]);
    ylim([0 5]);
end
xlabel('position (kb)');
ylabel('plectoneme size (kb)');
xlim([0 DNAinfo.DNAlen_bp/1000]);   % unit:kb
ylim([0 7]);
grid on;
hold off;

% in 2D pseudocolor
figure(752);clf
colormap('jet');
pt_connected=[];
for pti=1:N_pt
    pt_connected=[pt_connected; [pt.pos_in_bp{pti}; pt.size_in_bp{pti}]'];
end
size_vector=(0:0.1:10);
pos_size_2Dhist=hist3(pt_connected,{base_vector,size_vector});
imagesc(base_vector,size_vector,pos_size_2Dhist');axis xy;
xlabel('position (kb)');
ylabel('plectoneme size (kb)');
xlim([0 DNAinfo.DNAlen_bp/1000]);   % unit:kb
ylim([0 7]);
grid on;
hold off;


%% plectoneme density and occupancy in base position
figure(75);
subplot(2,1,2);
% build pt den
% pt_prof.density_weight_by_size=base_vector*0;
% pt_prof.density_no_weight=base_vector*0;
% 
% for pti=1:N_pt
%     c_pos=floor(1000*pt.pos_in_bp{pti}/base_resolution)+1;  % +1 as matlab vector starts from 1, unit kb
%     c_size=pt.size_in_bp{pti};
%     for tfi=1:length(c_pos)
%         pt_prof.density_weight_by_size(c_pos(tfi))=pt_prof.density_weight_by_size(c_pos(tfi))+c_size(tfi);
%         pt_prof.density_no_weight(c_pos(tfi))=pt_prof.density_no_weight(c_pos(tfi))+1;
%     end
% end
% 
% pt_prof.density_weight_by_size=pt_prof.density_weight_by_size/sum(pt_prof.density_weight_by_size(:));
% pt_prof.density_no_weight=pt_prof.density_no_weight/sum(pt_prof.density_no_weight(:));
% pt_prof.density_no_weight_raw=pt_prof.density_no_weight;

% plot(base_vector,pt_prof.density_weight_by_size,'r');hold on;
% plot(base_vector,pt_prof.density_no_weight,'b');hold off;
% [AX,H1,H2]=plotyy(base_vector_2x,pt_prof_2x.density_no_weight,base_vector_2x,pt_prof_2x.occupancy_inc_size);
[AX,H1,H2]=plotyy(base_vector,pt_prof.density_no_weight,base_vector,pt_prof.occupancy_inc_size);
legend('prob. (tip only)','ocup. (inc. size)','Location','best');

xlabel('position (kb)');
set(AX(1),'XLim',[0 DNAinfo.DNAlen_bp/1000],'YLim',[0 0.1])
set(AX(2),'XLim',[0 DNAinfo.DNAlen_bp/1000],'YLim',[0 1])
set(AX(2),'XTick',[],'YColor',[1 0 0]);
set(get(AX(1),'Ylabel'),'String','Probability density');
set(get(AX(2),'Ylabel'),'String','Occupancy (frame^-^1)');
set(H1,'Color',[0 0 0],'LineWidth',3);
set(H2,'Color',[1 0 0],'LineWidth',3);
grid on;
hold off;



%% calc nucleation rate
[pt.nuc_rate_avr, pt.nuc_rate]=update_nuc_rates(pt,initval,base_vector,Kymo_pt_occupation);

figure(721)
subplot(2,1,1);
% plot(pt.nuc_rate(:,1), pt.nuc_rate(:,2),'r');hold on;
% plot(pt.nuc_rate(:,1), pt.nuc_rate(:,3),'k');
plot(pt.nuc_rate(:,1), pt.nuc_rate(:,4),'g');hold on;
plot(pt.nuc_rate(:,1), pt.nuc_rate(:,5),'b');hold off;
xlim([0 DNAinfo.DNAlen_bp/1000]);
ylabel('Nucleation rate (s^-^1)');
xlabel('Position (kbp)');
legend('Nucleation','Termination');legend('boxoff');



%% calc flux

moved_L=[];
moved_R=[];
moved_S=[];
for pti=1:N_pt
    c_pt=pt.pos_in_bp{pti};
    c_pt_dspl=c_pt(1:end-1)-c_pt(2:end);
    
    moved_L=[moved_L c_pt(c_pt_dspl>0.5)];
    moved_R=[moved_R c_pt(c_pt_dspl<-0.5)];
    moved_S=[moved_S c_pt(-0.5 <= c_pt_dspl & c_pt_dspl <= 0.5)];
end

figure(721)
subplot(2,1,2);
hist_A=hist([moved_R moved_L moved_S],base_vector);
hist_L=hist(moved_L,base_vector)./hist_A;
hist_R=hist(moved_R,base_vector)./hist_A;
hist_S=hist(moved_S,base_vector)./hist_A;


plot(base_vector,hist_L,'r');hold on;
plot(base_vector,hist_S,'k');
plot(base_vector,hist_R,'b');hold off;
xlim([0 DNAinfo.DNAlen_bp/1000]);
xlabel('Position (kbp)');
ylabel('Probability');
legend('Left','Stay','Right','Location','Best');legend('boxoff');







%% plectoneme size vs. Diffusion coefficient
%
% figure(77);
% avr_pt_size=zeros(N_pt,1);
% avr_pt_position=zeros(N_pt,1);
% avr_pt_movement=zeros(N_pt,1);
% for pti=1:N_pt
%     c_pos=pt.pos_in_bp{pti};
%     c_size=pt.size_in_bp{pti};
%
%     avr_pt_size(pti)=mean(c_size);
%     avr_pt_position(pti)=mean(c_pos);
%     avr_pt_movement(pti)=std(c_pos);
% end
%
%
% % Size vs. Diff. coeff.
% subplot(2,2,1);
% plot(0);
% if n_valid_MSD~=0
%     plot(selected_size_mean,diffcoef_px,'.')
%     ylim([0 max(diffcoef_px)*1.05]);
% end
% xlim([0 3]);
% xlabel('Size (kb)','FontSize',12);
% ylabel('Diff. Coeff. (um^2/s)','FontSize',12);
% grid on;
%
% subplot(2,2,3);
% plot(avr_pt_size,avr_pt_movement,'.')
% xlabel('Size (bp)','FontSize',12);
% ylabel('STD of position (bp)','FontSize',12);
% grid on;
%
%
% % plectoneme lifetime vs. Diffusion coefficient
% subplot(2,2,2);
% plot(pt.pt_lifetime_all,avr_pt_size,'.')
% xlabel('Lifetime (s)','FontSize',12);
% ylabel('Size (bp)','FontSize',12);
% grid on;
% xlim([0 2]);
%
% subplot(2,2,4);
% plot(pt.pt_lifetime_all,avr_pt_movement,'.')
% xlabel('Lifetime (s)','FontSize',12);
% ylabel('STD of position (bp)','FontSize',12);
% grid on;
% xlim([0 2]);



% Save everything
save([kymo_output_dir '\allvariables3']);


