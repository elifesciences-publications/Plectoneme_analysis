function mxi=MG_Get1Dpeaks(h_profile,base_per_px,sigma,Npx2sum_sizedetermin)
%Find peaks in a vector h_profile

Npx=length(h_profile);

% find local max with window 5px
wdw=2;
all=[];
for pxi=wdw+1:Npx-wdw
    [~,mid]=max(h_profile(pxi-wdw:pxi+wdw));
    if mid==wdw+1
        all=[all pxi];
    end
end

% smooth h_profile
window_size=Npx2sum_sizedetermin;
new_h_profile=h_profile*0;
for pxi=window_size+1:Npx-window_size
    new_h_profile(pxi)=sum(h_profile(pxi-window_size:pxi+window_size)-base_per_px);
end

% select the peaks above the threshold
hi=new_h_profile(all)>base_per_px*sigma;
mxi=(all(hi));


% figure(511);
% plot(rij);hold on;
% plot(1:length(rij),rij*0+mn+sig*st);
% plot(mxi,rij(mxi),'k.');hold off;
% input('dksjdf');

end