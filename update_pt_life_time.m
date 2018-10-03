
function pt=update_pt_life_time(pt,initval,base_vector)

%% build individual plectoneme lifetime
pt.pt_lifetime_all=zeros(pt.N_pt,1);
for ttpti=1:pt.N_pt
    pt.pt_lifetime_all(ttpti)=length(pt.pos_in_px{ttpti})*initval.SecondsPerFrame;
end
pt.mean_pt_lifetime=mean(pt.pt_lifetime_all);

%% build position dependent lifetime

lifetimes_at_nuc_position=cell(length(base_vector),1);
lifetimes_at_ter_position=cell(length(base_vector),1);
lifetimes_at_mean_position=cell(length(base_vector),1);
base_vector2=base_vector-0.25;

for ttpti=1:pt.N_pt
    c_nuc_pos=pt.pos_in_bp{ttpti}(1);
    c_ter_pos=pt.pos_in_bp{ttpti}(end);
    c_mean_pos=mean(pt.pos_in_bp{ttpti});
    
    c_pt_lifetime=length(pt.pos_in_bp{ttpti})*initval.SecondsPerFrame;
    
    c_nuc_pos_id=find(base_vector2>c_nuc_pos,1,'first')-1;
    c_ter_pos_id=find(base_vector2>c_ter_pos,1,'first')-1;
    c_mean_pos_id=find(base_vector2>c_mean_pos,1,'first')-1;
    
    lifetimes_at_nuc_position{c_nuc_pos_id}=[lifetimes_at_nuc_position{c_nuc_pos_id} c_pt_lifetime];
    lifetimes_at_ter_position{c_ter_pos_id}=[lifetimes_at_ter_position{c_ter_pos_id} c_pt_lifetime];
    lifetimes_at_mean_position{c_mean_pos_id}=[lifetimes_at_mean_position{c_mean_pos_id} c_pt_lifetime];
end


pt.pt_lifetimes_at_nuc_position=base_vector*0;
pt.pt_lifetimes_at_ter_position=base_vector*0;
pt.pt_lifetimes_at_mean_position=base_vector*0;

for bvi=1:length(base_vector)
    pt.pt_lifetimes_at_nuc_position(bvi)=mean(lifetimes_at_nuc_position{bvi}(:));
    pt.pt_lifetimes_at_ter_position(bvi)=mean(lifetimes_at_ter_position{bvi}(:));
    pt.pt_lifetimes_at_mean_position(bvi)=mean(lifetimes_at_mean_position{bvi}(:));
end
