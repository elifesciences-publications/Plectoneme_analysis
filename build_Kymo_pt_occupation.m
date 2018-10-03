function [Kymo_pt_occupation, DNA_bins]=build_Kymo_pt_occupation(pt,DNAinfo,base_vector)

% unit_length_DNA=0.5;    % unit: kbp

% N_DNA_bin=floor(DNAinfo.DNAlen_bp/(unit_length_DNA*1000));
% DNA_bins=(1:N_DNA_bin)*unit_length_DNA-unit_length_DNA/2;

DNA_bins=base_vector;
N_DNA_bin=length(base_vector);
unit_length_DNA=DNA_bins(2)-DNA_bins(1);

Kymo_pt_occupation=zeros(N_DNA_bin,pt.N_frame_of_kymo_analyzed);

for ttpti=1:pt.N_pt
    for pt_tti=1:length(pt.pos_in_bp{ttpti})
        c_range1=floor( (pt.pos_in_bp{ttpti}(pt_tti)-pt.size_in_bp{ttpti}(pt_tti)/2)/unit_length_DNA )+1;
        c_range2=floor( (pt.pos_in_bp{ttpti}(pt_tti)+pt.size_in_bp{ttpti}(pt_tti)/2)/unit_length_DNA )+1;
        if c_range1 < 1
            c_range1 = 1;
        end
        if c_range2 > N_DNA_bin
            c_range2 = N_DNA_bin;
        end
        Kymo_pt_occupation(c_range1:c_range2,pt.start_fr(ttpti)+pt_tti-1)=1;
    end
end
% imagesc(base_vector-unit_length_DNA/2,1:pt.N_frame_of_kymo_analyzed,Kymo_pt_occupation');