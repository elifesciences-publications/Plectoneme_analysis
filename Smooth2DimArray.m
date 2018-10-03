
function kymo_smoothed=Smooth2DimArray(kymo_org,DNAinfo)
%Blur along X-direction
if nargin==1
    %% original code by Marijn => results in half sized kymo in x direction
    % kymo_smoothed=0;
    % for r=1:size(kymo_org,1)
    %     for k=1:length(kymo_org(r,:))
    %         if k==1
    %             kymo_smoothed(r,k)=mean(kymo_org(r,1:2));
    %         elseif k==length(kymo_org(r,:))
    %             kymo_smoothed(r,k)=mean(kymo_org(r,k-1:k));
    %         else
    %             kymo_smoothed(r,k)=(sum(kymo_org(r,k-1:k+1))+kymo_org(r,k))/4;     %% moving average filter with width 3, that counts the point itself twice
    %         end
    %     end
    % end

    %% Modifycation of SHKIM => Keep the x pixel number

    [Nframe,Npx]=size(kymo_org);
    kymo_smoothed=zeros(Nframe,Npx);

    % kymo_smoothed2=zeros(Nframe,Npx);
    for fri=1:Nframe
    %     kymo_smoothed(fri,:) = smooth(kymo_org(fri,:),0.1,'rloess');  % the method used is efficient to remove outliers
        kymo_smoothed(fri,:) = smooth(kymo_org(fri,:),3,'moving');
    %     figure(531)
    %     plot(kymo_org(fri,:),'o');   hold on;
    %     plot(kymo_smoothed(fri,:));  hold off;
    %     ylim([-6 10]);
    %     title(fri);
    %     input('whatever');
    end
elseif nargin==2
    disp('tough');
    disp('the problem will be fixed later...');
end