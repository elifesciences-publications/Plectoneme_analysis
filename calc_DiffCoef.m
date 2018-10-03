
function [pt,msd,dsp]=calc_DiffCoef(pt,initval,min_trace_len_ana,min_pt_size_ana)
%% input arguments
% min_trace_len=10;   % discard diffusion shorter than min_trace_len
% min_pt_size_ana=0.3;  % unit: kb

%% calc displacement
for pt_i=1:pt.N_pt
    dsp.px(pt_i)={(pt.pos_in_px{pt_i}(:)-pt.pos_in_px{pt_i}(1))*initval.Px2um};
    dsp.bp(pt_i)={(pt.pos_in_bp{pt_i}(:)-pt.pos_in_bp{pt_i}(1))};
    dsp.times(pt_i)={(1:length(dsp.px{pt_i}))*initval.SecondsPerFrame};
end

%% Calc MSD
for pt_i=1:pt.N_pt
    c_pt_um=pt.pos_in_px{pt_i}*initval.Px2um;
    c_pt_bp=pt.pos_in_bp{pt_i};
    c_msd_um=nan(1,length(c_pt_um)-1);
    c_msd_bp=nan(1,length(c_pt_um)-1);
    for tau=1:length(c_pt_um)-1
        c_msd_um(tau)=mean( (c_pt_um(1:end-tau)-c_pt_um(1+tau:end)).^2 );
        c_msd_bp(tau)=mean( (c_pt_bp(1:end-tau)-c_pt_bp(1+tau:end)).^2 );
    end
    msd.px(pt_i)={c_msd_um};
    msd.bp(pt_i)={c_msd_bp};
    msd.times(pt_i)={(1:length(msd.px{pt_i}))*initval.SecondsPerFrame};
end

%% Convert msd and displacement to matrics form
clear('msd_easy_px','msd_easy_bp','dsp.px_easy','dsp.bp_easy','msd_times_easy','dsp_times_easy');
for pt_i=1:pt.N_pt
    msd.easy_px(1:length(msd.px{pt_i}),pt_i)=msd.px{pt_i};
    msd.easy_bp(1:length(msd.bp{pt_i}),pt_i)=msd.bp{pt_i};

    dsp.px_easy(1:length(dsp.px{pt_i}),pt_i)=dsp.px{pt_i};
    dsp.bp_easy(1:length(dsp.bp{pt_i}),pt_i)=dsp.bp{pt_i};

    msd.times_easy(1:length(msd.px{pt_i}),pt_i)=msd.times{pt_i};
    dsp.times_easy(1:length(dsp.px{pt_i}),pt_i)=dsp.times{pt_i};
end
%%
valid_id=msd.times_easy==0;
msd.times_easy(valid_id)=NaN;
msd.easy_px(valid_id)=NaN;
msd.easy_bp(valid_id)=NaN;
dsp.px_easy(dsp.times_easy==0)=NaN;
dsp.bp_easy(dsp.times_easy==0)=NaN;
dsp.times_easy(dsp.times_easy==0)=NaN;



%% calc diffusion coefficient by linear fit of data


[tmax,~]=size(msd.times_easy);
time=(1:tmax)'*initval.SecondsPerFrame;
n_valid_MSD=0;
msd.analyzed_px=[];
msd.analyzed_bp=[];
msd.analyzed_time=[];
diffcoef_px=[];
diffcoef_bp=[];
selected_size_mean=[];

for pt_i=1:pt.N_pt
    %% pick one plectonme trace
%     c_msd_px=msd_px{pt_i};
%     c_msd_bp=c_msd_bp{pt_i};
    c_msd_um=msd.easy_px(:,pt_i);
    c_msd_bp=msd.easy_bp(:,pt_i);
    c_size=mean(pt.size_in_bp{pt_i});
    
    % remove trailing nan
    valid_id_px=~isnan(c_msd_um);
    c_msd_um=c_msd_um(valid_id_px);
    c_msd_bp=c_msd_bp(valid_id_px);
    
    c_time=time(valid_id_px);
    c_time_max(pt_i)=c_time(end);
    
    if sum(valid_id_px) > min_trace_len_ana && c_size > min_pt_size_ana
        n_valid_MSD=n_valid_MSD+1;

        selected_size_mean(n_valid_MSD)=c_size;
        
        [tot_N,~]=size(c_msd_um);
        W_vecter=(1:tot_N);
        
        %% fit the MSD in um
        a0=(c_msd_um(end)-c_msd_um(1))/(c_time(end)-c_time(1));
        
        s2 = fitoptions('Method','NonlinearLeastSquares',...
            'Lower', 0,'Upper',Inf,...
            'Startpoint',a0);
        f2 = fittype('a*x','options',s2);
        
        sigma_msd = 2/3*W_vecter .* ( tot_N+1 - W_vecter );
        Weights=1./sigma_msd;   % Weight = 1/std^2 or 1/sigma
        
        [c2_px,gof2] = fit(c_time,c_msd_um,f2,'Weight',Weights);
        diffcoef_px(n_valid_MSD)=c2_px.a;
        
        msd.analyzed_px((1:length(c_msd_um)),n_valid_MSD)=c_msd_um;
        msd.analyzed_time((1:length(c_msd_um)),n_valid_MSD)=(1:length(c_msd_um))*initval.SecondsPerFrame;
        
        %% fit the MSD in bp
        a0=(c_msd_bp(end)-c_msd_bp(1))/(c_time(end)-c_time(1));
        
        s2 = fitoptions('Method','NonlinearLeastSquares',...
            'Lower', 0,'Upper',Inf,...
            'Startpoint',a0);

        [c2_bp,gof2] = fit(c_time,c_msd_bp,f2,'Weight',Weights);
        diffcoef_bp(n_valid_MSD)=c2_bp.a;
        
        msd.analyzed_bp((1:length(c_msd_bp)),n_valid_MSD)=c_msd_bp;
        
    else
        
    end
end
%%


if ~isempty(msd.analyzed_px)
    % put the result into pt
	pt.D_mean_px=mean(diffcoef_px);
    pt.D_mean_bp=mean(diffcoef_bp);
    pt.D_std_px=std(diffcoef_px);
    pt.D_std_bp=std(diffcoef_bp);
    pt.D_se_px=pt.D_std_px/sqrt(n_valid_MSD);
    pt.D_se_bp=pt.D_std_bp/sqrt(n_valid_MSD);

    % calculate averaged MSD
    msd.analyzed_px(msd.analyzed_time==0)=NaN;
    msd.analyzed_bp(msd.analyzed_time==0)=NaN;
    msd.analyzed_time(msd.analyzed_time==0)=NaN;
    msd.mean_msd_px=mean(msd.analyzed_px,2);
    msd.mean_msd_bp=mean(msd.analyzed_bp,2);
    
    valid_id=~isnan(msd.mean_msd_px);
    msd.mean_msd_px=msd.mean_msd_px(valid_id);
    msd.mean_msd_bp=msd.mean_msd_bp(valid_id);
    
    msd.analyzed_se_px=std(msd.analyzed_px,0,2)/sqrt(n_valid_MSD);
    msd.analyzed_se_bp=std(msd.analyzed_bp,0,2)/sqrt(n_valid_MSD);
    msd.analyzed_se_px=msd.analyzed_se_px(valid_id);
    msd.analyzed_se_bp=msd.analyzed_se_bp(valid_id);
    
    % fit the averaged MSD in px
    c_time=(1:length(msd.mean_msd_px))'*initval.SecondsPerFrame;
    a0=(msd.mean_msd_px(end)-msd.mean_msd_px(1))/(c_time(end)-c_time(1));
    s2 = fitoptions('Method','NonlinearLeastSquares',...
        'Lower', 0,'Upper',Inf,...
        'Startpoint',a0);
    s3 = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0],'Upper',[Inf,Inf],...
        'Startpoint',[a0,0.9]);
    f2 = fittype('a*x','options',s2);
    f3 = fittype('a*x^n','options',s3);
    [tot_N,~]=size(msd.mean_msd_px);
    Weights=msd.analyzed_se_px;
    
    [c2_px,gof2] = fit(c_time,msd.mean_msd_px,f2,'Weight',Weights);
    [c3_px,gof3] = fit(c_time,msd.mean_msd_px,f3,'Weight',Weights);

    pt.D_of_meanMSD_freeDiff_px=c2_px.a;
    pt.D_of_meanMSD_SubDiff_px=c3_px.a;
    pt.n_of_meanMSD_SubDiff_px=c3_px.n;


    % fit the averaged MSD in bp
    a0=(msd.mean_msd_bp(end)-msd.mean_msd_bp(1))/(c_time(end)-c_time(1));
    s2 = fitoptions('Method','NonlinearLeastSquares',...
        'Lower', 0,'Upper',Inf,...
        'Startpoint',a0);
    s3 = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0],'Upper',[Inf,Inf],...
        'Startpoint',[a0,0.9]);
    
    Weights=msd.analyzed_se_bp;
    
    [c2_bp,gof2] = fit(c_time,msd.mean_msd_bp,f2,'Weight',Weights);
    [c3_bp,gof3] = fit(c_time,msd.mean_msd_bp,f3,'Weight',Weights);

    pt.D_of_meanMSD_freeDiff_bp=c2_bp.a;
    pt.D_of_meanMSD_SubDiff_bp=c3_bp.a;
    pt.n_of_meanMSD_SubDiff_bp=c3_bp.n;
    
    msd.msd_fit_linear=c2_bp;
    msd.msd_fit_sub_diffusion=c3_bp;
    msd.msd_fit_linear_px=c2_px;
    msd.msd_fit_sub_diffusion_px=c3_px;
    msd.n_valid_MSD=n_valid_MSD;
    msd.tmax=tmax;
    msd.min_trace_len_ana=min_trace_len_ana;
    msd.min_pt_size_ana=min_pt_size_ana;
end



















