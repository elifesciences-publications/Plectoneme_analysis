
function [pt_prof]=build_ptden(pt,base_vector,Kymo_pt_occupation)


%% calc pt density (probability)
pt_prof.density_weight_by_size=base_vector*0;
pt_prof.hist_tip=base_vector*0;

base_resolution=base_vector(2)-base_vector(1);

for pti=1:pt.N_pt
    c_pos=round(pt.pos_in_bp{pti}/base_resolution)+1;  % +1 as matlab vector starts from 1, unit kb
    c_size=pt.size_in_bp{pti};
    for tfi=1:length(c_pos)
        pt_prof.density_weight_by_size(c_pos(tfi))=pt_prof.density_weight_by_size(c_pos(tfi))+c_size(tfi);
        pt_prof.hist_tip(c_pos(tfi))=pt_prof.hist_tip(c_pos(tfi))+1;
    end
end

pt_prof.density_weight_by_size=pt_prof.density_weight_by_size/sum(pt_prof.density_weight_by_size(:));
pt_prof.density_no_weight=pt_prof.hist_tip/sum(pt_prof.hist_tip(:));


%% calc occupancy
pt_prof.hist_ptmicDNA=sum(Kymo_pt_occupation,2);

pt_prof.occupancy_exc_size=pt_prof.hist_tip/pt.N_frame_of_kymo_analyzed;
pt_prof.occupancy_inc_size=pt_prof.hist_ptmicDNA/pt.N_frame_of_kymo_analyzed;






