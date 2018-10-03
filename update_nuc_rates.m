function [nuc_rate_avr, nuc_rates]=update_nuc_rates(pt,initval,base_vector,Kymo_pt_occupation)
%% average nuc_rate of the kymograph
nuc_rate_avr=(pt.N_pt-sum(pt.start_fr==1))/(pt.N_frame_of_kymo_analyzed*initval.SecondsPerFrame);  % Unit: per second


%% position dependent nuc_rate method1
nuc_rates(:,1)=base_vector;

for pti=1:pt.N_pt
    c_nuc_pos(pti,1)=pt.pos_in_bp{pti}(1);
    c_nuc_pos(pti,2)=pt.pos_in_bp{pti}(end);
end

nuc_rate_raw(:,1)=hist(c_nuc_pos(:,1),nuc_rates(:,1));
nuc_rate_raw(:,2)=hist(c_nuc_pos(:,2),nuc_rates(:,1));

% nuc_rates(:,2)=nuc_rate_raw(:,1)./(initval.SecondsPerFrame*(pt.N_frame_of_kymo_analyzed - pt_den_in_base_no_weight_raw + nuc_rate_raw(:,1)));
% nuc_rates(:,3)=nuc_rate_raw(:,2)./(initval.SecondsPerFrame*(pt.N_frame_of_kymo_analyzed - pt_den_in_base_no_weight_raw + nuc_rate_raw(:,2)));
nuc_rates(:,2)=nuc_rate_raw(:,1);
nuc_rates(:,3)=nuc_rate_raw(:,2);

%% position dependent nuc_rate method2

available_frames=sum(~Kymo_pt_occupation,2);
occupied_frames=sum(Kymo_pt_occupation,2);

nuc_rates(:,4)=nuc_rate_raw(:,1)./(initval.SecondsPerFrame*(available_frames + nuc_rate_raw(:,1)));
% nuc_rates(:,5)=nuc_rate_raw(:,2)./(initval.SecondsPerFrame*(available_frames + nuc_rate_raw(:,2)));
nuc_rates(:,5)=nuc_rate_raw(:,2)./(initval.SecondsPerFrame*(occupied_frames));

