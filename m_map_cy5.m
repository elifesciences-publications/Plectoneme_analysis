function ImRGB_mapped=m_map_cy5(ImCy3,ImCy5,initval)
% function m_map_cy5(ImCy3,ImCy5,initval)
if nargin<3 % development purpose
    %     ImCy3=[];
    %     ImCy5=[];
    debug_mode=1;
else
    debug_mode=0;
end

%% load mapping file
% check if map is already maden
if ~debug_mode
    path_for_map=initval.plecdir;
    token_id=strfind(path_for_map,'\');
    path_for_map=path_for_map(1:token_id(end));
    
    use_oldmap=1;
    if exist([path_for_map 'map.mapcoef'],'file')
        %     user_ans=menu(['Map found in [' path_for_map(1:end-1) '].'],'Load old','Make new');
        disp(' ');
        disp(['Map found in [' path_for_map(1:end-1) '].']);
        user_ans=input('1. Load old, 2. Make new [default=1]  ');
        if user_ans==2
            use_oldmap=2;
        end
    else
        use_oldmap=2;
    end
    
    if use_oldmap == 1    % load existing map
        map_coef=load_map(path_for_map);
    else    %% make new mapping file
        [map_coef,~,~]=mapcreator(path_for_map);
    end
elseif debug_mode
    [map_coef,ImCy3,ImCy5]=mapcreator(cd());
end

%% make mapped image
ImRGB_mapped=mergeimage_RGB(ImCy3, ImCy5, map_coef);

end


function [map,ImGreen,ImRed]=mapcreator(path_for_map_save)

%% get bead images
dir_green_bead=uigetdir('K:\bn\alg\Shared\Mahipal-Kim\Beads for mapping','Green beads');
dir_red_bead=uigetdir('K:\bn\alg\Shared\Mahipal-Kim\Beads for mapping','Red beads');

green_filelist = dir([dir_green_bead '\*.tiff']);
[Nimages_green,~]=size(green_filelist);

red_filelist = dir([dir_red_bead '\*.tiff']);
[Nimages_red,~]=size(red_filelist);

% read green bead images
for fr_id=1:Nimages_green
    TifLink = Tiff([dir_green_bead '\' green_filelist(fr_id).name] , 'r');
    Im=double(TifLink.read());
    Im=Im-median(Im(:));
    Im=Im./max(Im(:));
    if fr_id==1
        ImGreen=Im;
    else
        ImGreen=ImGreen+Im;
    end
end
ImGreen=ImGreen/Nimages_green;

% read red bead images
for fr_id=1:Nimages_red
    TifLink = Tiff([dir_red_bead '\' red_filelist(fr_id).name] , 'r');
    Im=double(TifLink.read());
    Im=Im-median(Im(:));
    Im=Im./max(Im(:));
    if fr_id==1
        ImRed=Im;
    else
        ImRed=ImRed+Im;
    end
end
ImRed=ImRed/Nimages_red;

%% Parameters
[XpixelSize,YpixelSize]=size(ImGreen);
framesize=XpixelSize*YpixelSize;
% NumFrameProc=10;
use_gauss_filter= 1; % 0:off, 1: on
circle_color=.5;
localarea_size=5;

frame(1:YpixelSize,1:XpixelSize)=ImGreen;
frame(1:YpixelSize,XpixelSize+1:XpixelSize*2)=ImRed;

%% get position of 5 molecules (10pt) for mapping
circled_frame=frame;
circled_frame(:,XpixelSize:XpixelSize+1)=1;   % add a line bet. donor and acceptor channel

figure('Units','pixels','Position',[10 50 XpixelSize*2 YpixelSize]);    hdl_mcfig=gcf;
axes('position',[0 0 1 1]);
imagesc(circled_frame);caxis([0 .8]);
axis image;
colormap('jet');

flg_pt=zeros(10,2);  % odd numbered raw for donor, even for acceptor
for i=1:10  % get positions of five molecules
    [flg_pt(i,2), flg_pt(i,1)]=ginput(1);
    flg_pt=floor(flg_pt);
    
    % find local maxima
    subframe=frame(flg_pt(i,1)-localarea_size:flg_pt(i,1)+localarea_size,...
        flg_pt(i,2)-localarea_size:flg_pt(i,2)+localarea_size);
    [tmp, tmpindex1]=max(subframe);
    [~, tmpindex2]=max(tmp);
    flg_pt(i,1)=tmpindex1(tmpindex2)-localarea_size-1+flg_pt(i,1);
    flg_pt(i,2)=tmpindex2-localarea_size-1+flg_pt(i,2);
    
    % draw circle to the selected molecule
    circled_frame=add_circle(circled_frame,flg_pt(i,1),flg_pt(i,2),circle_color);
    imagesc(circled_frame);caxis([0 .8]);    axis image;drawnow;
    
    % Half-frame problem correction
    if mod(i,2)==0
        flg_pt(i,2)=flg_pt(i,2)-XpixelSize;
    end
end

%% Get temporary map coef. with the 3 molecules
fit_poly_order=2;   % this will be used below...
x=[flg_pt(1,:) flg_pt(3,:) flg_pt(5,:) flg_pt(7,:) flg_pt(9,:)]';
y=[flg_pt(2,:) flg_pt(4,:) flg_pt(6,:) flg_pt(8,:) flg_pt(10,:)]';
ft = fittype( 'mapfn3mol_2nd( x, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4)' );
f = fit( x, y, ft, 'StartPoint', [0, 0, 1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1] );

%% find molecules
[mol_pos_G, num_mol_G]=peakfinderembeded(ImGreen,localarea_size,1,use_gauss_filter);
[mol_pos_R, num_mol_R]=peakfinderembeded(ImRed,localarea_size,1,use_gauss_filter);
disp([' green beads: ' num2str(num_mol_G) ', red beads: ' num2str(num_mol_R)]);

%% Do multiple round of poly-fit
num_paired_prev=-1;
for ii=1:25
    %% Find the pair of all the molecules with temporary map coef.
    [pairsfound,num_paired,mapped_frame]=find_pairs(frame,mol_pos_G,num_mol_G,mol_pos_R,num_mol_R,f,fit_poly_order,circle_color,XpixelSize);
    disp([num2str(num_paired) ' pairs found in ' num2str(ii) 'th iteration']);
    if num_paired > 0.9*num_mol_G
        disp('Fit statisfied');
        break;
    elseif num_paired == num_paired_prev
        disp('Fit accepted.');
        break;  
    end
    num_paired_prev=num_paired;
    %% Get map coef with the pairs found in 3-moleule-fit
    for i=1:num_paired
        x(i*2-1,:)=pairsfound(i,1);
        x(i*2,:)=pairsfound(i,2);
        y(i*2-1,:)=pairsfound(i,3);
        y(i*2,:)=pairsfound(i,4);
    end
    fit_poly_order=4;   % n-th order polynomial
    f=fit_poly(x,y,fit_poly_order);
end


%% save mapping coefficient in dir_green_bead
clear('map')
writefilename=[path_for_map_save 'map.mapcoef'];
mapfid=fopen(writefilename,'w');
A=coeffnames(f);
B=coeffvalues(f);
for i=1:size(A);
    fprintf(mapfid,'%s %g\n',A{i},B(i));
end
map(1,:)=B(1:size(A)/2);
map(2,:)=B(size(A)/2+1:size(A));
fclose(mapfid);

%% plot image file
figure(hdl_mcfig);
mapped_frame(:,XpixelSize:XpixelSize+1)=1;   % add a line bet. donor and acceptor channel
imagesc(mapped_frame);caxis([0 0.8]);
figure();subplot('Position',[0 0 1 1]);
ImRGB_mapped=mergeimage_RGB(ImGreen, ImRed, map);
imshow(ImRGB_mapped*3);

% pause(2);
% close(hdl_mcfig);
end


function f=fit_poly(x,y,fit_poly_order)
% fit_poly_order=4;   % n-th order polynomial
%
if fit_poly_order==3    % at least 7 molecules required
    ft = fittype( 'mapfn( x, a0, a1, a2, a3, a4, a5, a6, b0, b1, b2, b3, b4, b5, b6)' );
    f = fit( x, y, ft, 'StartPoint', [0, 1, 0.1, 0.1, 0.1 0.01, 0.01, 0, 1, 0.1, 0.1, 0.1 0.01, 0.01] );
elseif fit_poly_order==4    % at least 9 molecules required
    ft = fittype( 'mapfn_4th( x, a0, a1, a2, a3, a4, a5, a6, a7, a8, b0, b1, b2, b3, b4, b5, b6, b7, b8)' );
    f = fit( x, y, ft, 'StartPoint', [0, 1, 0.1, 0.1, 0.1 0.01, 0.01, 0.01, 0.01, 0.1, 1, 0.1, 0.1, 0.1 0.01, 0.01, 0.01, 0.01] );
end
end

function [pairsfound,num_paired,mapped_frame]=find_pairs(frame,mol_pos_G,num_mol_G,mol_pos_R,num_mol_R,f,fit_poly_order,circle_color,XpixelSize)
num_paired=0;
pairsfound=zeros(1,4);
mapped_frame=frame;

for i=1:num_mol_G
    if fit_poly_order==2
        % 2nd order (at least 5 molecule required)
        mappedx= round(f.a0 + f.a1*mol_pos_G(i,1) + f.a2*mol_pos_G(i,2) + f.a3*mol_pos_G(i,1)^2 + f.a4*mol_pos_G(i,2)^2);
        mappedy= round(f.b0 + f.b1*mol_pos_G(i,1) + f.b2*mol_pos_G(i,2) + f.b3*mol_pos_G(i,1)^2 + f.b4*mol_pos_G(i,2)^2);
    elseif fit_poly_order==3
        mappedx= round(f.a0 + f.a1*mol_pos_G(i,1) + f.a2*mol_pos_G(i,2) + f.a3*mol_pos_G(i,1)^2 + f.a4*mol_pos_G(i,2)^2 + f.a5*mol_pos_G(i,1)^3 + f.a6*mol_pos_G(i,2)^3);
        mappedy= round(f.b0 + f.b1*mol_pos_G(i,1) + f.b2*mol_pos_G(i,2) + f.b3*mol_pos_G(i,1)^2 + f.b4*mol_pos_G(i,2)^2 + f.b5*mol_pos_G(i,1)^3 + f.b6*mol_pos_G(i,2)^3);
    elseif fit_poly_order==4
        mappedx= round(f.a0 + f.a1*mol_pos_G(i,1) + f.a2*mol_pos_G(i,2) + f.a3*mol_pos_G(i,1)^2 + f.a4*mol_pos_G(i,2)^2 + f.a5*mol_pos_G(i,1)^3 + f.a6*mol_pos_G(i,2)^3 + f.a7*mol_pos_G(i,1)^4 + f.a8*mol_pos_G(i,2)^4);
        mappedy= round(f.b0 + f.b1*mol_pos_G(i,1) + f.b2*mol_pos_G(i,2) + f.b3*mol_pos_G(i,1)^2 + f.b4*mol_pos_G(i,2)^2 + f.b5*mol_pos_G(i,1)^3 + f.b6*mol_pos_G(i,2)^3 + f.b7*mol_pos_G(i,1)^4 + f.b8*mol_pos_G(i,2)^4);
    end
    
    % check if there is acceptor at the mapped position
    for k=1:num_mol_R
        %             if (mol_pos(k,1)==mappedx && mol_pos(k,2)==mappedy),
        if (mol_pos_R(k,1)>mappedx-3 && mol_pos_R(k,1)<mappedx+3 && mol_pos_R(k,2)>mappedy-3 && mol_pos_R(k,2)<mappedy+3)
            %                 disp(['(' num2str(mol_pos_G(i,1)) ',' num2str(mol_pos_G(i,2)) ')(' num2str(mappedx) ',' num2str(mappedy) ')(' num2str(mol_pos_R(k,1)) ',' num2str(mol_pos_R(k,2)) ')']);
            
            num_paired=num_paired+1;
            pairsfound(num_paired,1)=mol_pos_G(i,1);
            pairsfound(num_paired,2)=mol_pos_G(i,2);
            pairsfound(num_paired,3)=mol_pos_R(k,1);
            pairsfound(num_paired,4)=mol_pos_R(k,2);
            
            mapped_frame=add_circle(mapped_frame,mol_pos_G(i,1),mol_pos_G(i,2),circle_color);
            mapped_frame=add_circle(mapped_frame,mol_pos_R(k,1),mol_pos_R(k,2)+XpixelSize,circle_color);
            %                 mapped_frame=add_circle(mapped_frame,mappedx,mappedy,circle_color);
        end
    end
end
% pairsfound=pairsfound(1:num_paired,:);
% disp([num2str(num_paired) ' molecules used for mapping']);

end


function map=load_map(path_for_map)
mapfid=fopen([path_for_map 'map.mapcoef'],'r');

i=0;
while(~isempty(fscanf(mapfid,'%s',1)))
    i=i+1;
    tmp(i)=fscanf(mapfid,'%g',1);
end
map(1,:)=tmp(1:length(tmp)/2);
map(2,:)=tmp(length(tmp)/2+1:end);
fclose(mapfid);

end


function ImRGB_mapped=mergeimage_RGB(D_image, A_image, map)

[X_size, Y_size]=size(D_image);
ImRGB_mapped=zeros(X_size,Y_size);
[~,tmp]=size(map);
if tmp==5   % 2nd order
    fit_poly_order=2;
elseif tmp==7	% third order
    fit_poly_order=3;
elseif tmp==9   % forth order
    fit_poly_order=4;
%     disp('4th order used');
end


for x=1:X_size
    for y=1:Y_size
        if fit_poly_order==2
            mappedx= round(map(1,1) +...
                map(1,2)*x + map(1,3)*y +...
                map(1,4)*x^2 + map(1,5)*y^2);
            mappedy= round(map(2,1) +...
                map(2,2)*x + map(2,3)*y +...
                map(2,4)*x^2 + map(2,5)*y^2);
        elseif fit_poly_order==3
            mappedx= round(map(1,1) +...
                map(1,2)*x + map(1,3)*y +...
                map(1,4)*x^2 + map(1,5)*y^2 +...
                map(1,6)*x^3 + map(1,7)*y^3);
            mappedy= round(map(2,1) +...
                map(2,2)*x + map(2,3)*y +...
                map(2,4)*x^2 + map(2,5)*y^2 +...
                map(2,6)*x^3 + map(2,7)*y^3);
        elseif fit_poly_order==4
            mappedx= round(map(1,1) +...
                map(1,2)*x + map(1,3)*y +...
                map(1,4)*x^2 + map(1,5)*y^2 +...
                map(1,6)*x^3 + map(1,7)*y^3 +...
                map(1,8)*x^4 + map(1,9)*y^4);
            mappedy= round(map(2,1) +...
                map(2,2)*x + map(2,3)*y +...
                map(2,4)*x^2 + map(2,5)*y^2 +...
                map(2,6)*x^3 + map(2,7)*y^3 +...
                map(2,8)*x^4 + map(2,9)*y^4);
        end
        
        %         mappedy=mappedy-Y_size;
        if mappedx < 1, mappedx = 1; end
        if mappedy < 1, mappedy = 1; end
        if mappedx > X_size, mappedx = X_size; end
        if mappedy > Y_size, mappedy = Y_size; end
        %         Intensity_merged(x,y)=D_image(x,y)+A_image(mappedx,mappedy);
        %         FRET_merged(x,y)=A_image(mappedx,mappedy)/(D_image(x,y)+A_image(mappedx,mappedy));
        ImRGB_mapped(x,y,1)=A_image(mappedx,mappedy);
        ImRGB_mapped(x,y,2)=D_image(x,y);
        ImRGB_mapped(x,y,3)=0;
    end
end

end

function frame=add_circle(frame, i, j,circle_color)
% circle_color=170;
[Xsize, Ysize]=size(frame);
if isscalar(i)
    if i>3 && i < Xsize-3 && j>3 && j<Ysize-3
        frame(i-3,j-1)=circle_color;
        frame(i-3,j+1)=circle_color;
        frame(i-1,j-3)=circle_color;
        frame(i-1,j+3)=circle_color;
        frame(i+1,j-3)=circle_color;
        frame(i+1,j+3)=circle_color;
        frame(i+3,j-1)=circle_color;
        frame(i+3,j+1)=circle_color;
    end
else
    num_mol=size(i);
    iarray=i;
    jarray=j;
    for curr_mol=1:num_mol
        i=iarray(curr_mol);
        j=jarray(curr_mol);
        if i>3 && i < Xsize-3 && j>3 && j<Ysize-3
            frame(i-3,j-1)=circle_color;
            frame(i-3,j+1)=circle_color;
            frame(i-1,j-3)=circle_color;
            frame(i-1,j+3)=circle_color;
            frame(i+1,j-3)=circle_color;
            frame(i+1,j+3)=circle_color;
            frame(i+3,j-1)=circle_color;
            frame(i+3,j+1)=circle_color;
        end
    end
end
end

function [mol_pos, num_mol, background]=peakfinderembeded(frame,localarea_size,threshold_index,use_gauss_filter)
%% parameters
border=localarea_size+1;
box_size=localarea_size*2+1;
[Xsize, Ysize]=size(frame);
num_mol=0;

%% make border box
borderbox=zeros(box_size);
borderbox(1,:)=1;
borderbox(box_size,:)=1;
borderbox(:,1)=1;
borderbox(:,box_size)=1;
centerbox=ones(box_size)-borderbox;
% borderbox_area=box_size*4-4;
borderbox_ind=borderbox(:)==1;
centerbox_area=(box_size-2)^2;

%% Create 2D-Gaussian Filtering
if (use_gauss_filter)
    % parameters for 2D-gaussian filtering
    sigma=3;
    
    % Make Gaussian Mask
    FixedGmask=zeros(box_size,box_size);
    for i=1:box_size
        for j=1:box_size
            FixedGmask(j,i)=exp( -((i-localarea_size-1).^2+(j-localarea_size-1).^2) ./ sigma );
        end
    end
    GmaskAmp = sum(sum(FixedGmask));
    
    % Convolute frame with 2D-Guassian Mask
    frameConvo=zeros(Xsize,Ysize);
    for i=border:Xsize-border
        for j=border:Ysize-border
            subframe=frame(i-localarea_size:i+localarea_size,...
                j-localarea_size:j+localarea_size);
            frameConvo(i,j)=sum(sum(subframe.*FixedGmask))/GmaskAmp;
        end
    end
    frame=frameConvo;
end

%% Find Local Maxima
for i=border:Xsize-border
    for j=border:Ysize-border
        subframe=frame(i-localarea_size:i+localarea_size,j-localarea_size:j+localarea_size);
        cur_pixel=subframe(localarea_size+1,localarea_size+1);
        
        [max_valX max_index]=max(subframe);
        [~, max_indexY]=max(max_valX);
        
        if max_index(max_indexY)==localarea_size+1 && max_indexY==localarea_size+1,
            % is local maxima?
            cur_borderbox=subframe(borderbox_ind);
            cur_background=mean(cur_borderbox);
            cur_bordermax=max(cur_borderbox);
            cur_borderbox_std=std(cur_borderbox);
            
            cur_intensity=sum(sum(subframe.*centerbox))/centerbox_area;
            
            if cur_pixel-cur_background > (cur_bordermax-cur_background)*3,
                if cur_intensity-cur_background > cur_borderbox_std*threshold_index,
                    % accept as a molecule if the intensity of the centor box
                    % is larger than background*threshold_index
                    num_mol=num_mol+1;
                    mol_pos(num_mol,1)=i;
                    mol_pos(num_mol,2)=j;
                    background(num_mol)=cur_background;
                end
            end
        end
    end
end
end
