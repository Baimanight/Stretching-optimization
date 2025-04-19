%% ---Stretching---

% Based on MSNoise1.6 python code 
% Fold double sides as computed waveform 
% Gauss fitting center difference between Ref and Days --> offset
% offset cannot be directly considered as clock error
% 
% fix a small bug and improve the smooth 

%-------------------
warning('off'); 
% python warning   map_coordinates contains a extra filter  
% the filter's code seems under cover    
clear;
%%------Setup-------

path = '/home/baimanight/Test';

STACK = '/STACKS';
directory = '/Str_smooth';
filter = '/01'; 
stack = '/005_DAYS';

REF = strcat(path,STACK,filter,'/REF');
pathi = strcat(path,STACK,filter,stack);
patho = strcat(path,directory,filter,stack);
path_offset = strcat(path,directory,filter);;
% path_offset/data_offset
% Each line corresponds to a pair of stations.
% like, RR_A-B; RR_A-C; RR_B-C; TT_A-B; TT_A-C...

maxoffset = 20;    % max oscillating offset
maxdiff = 10;      % max difference of oscillating offsets
dura_lag = 7;      % min elimating dura
sample_rate = 20;
sidewidth = 120;
minlag = 15;
maxlag = 65;
st_dura = 0.02;    % dvv dura
step = 0.001;      % min step  
st_num = fix(1/step);

%----Setup over----
if ~exist(strcat(path,directory), 'dir')
    mkdir(strcat(path,directory));
end
if ~exist(strcat(path,directory,'/01'), 'dir')
    mkdir(strcat(path,directory,'/01'));
end
if ~exist(strcat(path,directory,'/01','/001_DAYS'), 'dir')
    mkdir(strcat(path,directory,'/01','/001_DAYS'));
end

offset_num = 1;
Ref = getdir(REF);
Stack = getdir(pathi);
lag_fullwidth = [-sidewidth,sidewidth];
Rightlag = (min(lag_fullwidth):1/sample_rate:max(lag_fullwidth))'; 

ndimage = py.importlib.import_module('scipy.ndimage');

mid = max(lag_fullwidth)*sample_rate;
replace_lag = zeros(1,2*mid+1);
crep1 = replace_lag; crep2 = replace_lag;
extra_num = fix((maxlag-minlag)*sample_rate*st_dura+10);  % for External interpolation

replace_lag(mid+minlag*sample_rate-extra_num:mid+maxlag*sample_rate+extra_num) = 1;
replace_lag(mid-maxlag*sample_rate-extra_num:mid-minlag*sample_rate+extra_num) = 1;
crep1(mid+minlag*sample_rate:mid+maxlag*sample_rate) = 1; crepr = crep1;
crep2(mid-maxlag*sample_rate:mid-minlag*sample_rate) = 1; crepl = crep2;


%% ----Loop---

for component_num = 1:length(Ref) 
    com_name = dir(pathi);comp_name = com_name(component_num+2).name;
    progressing = sprintf('%d/%d',component_num,length(Ref));
    display(strcat('Computing Component ',sprintf(' %s  %s',comp_name, ...
        progressing)));      disp(' ');
    Ref_pair = getdir(Ref{component_num});
    Right_pair = getdir(Stack{component_num});
    pai_name=dir(Stack{component_num});

    if ~exist(strcat(patho,'/',comp_name), 'dir')
     mkdir(strcat(patho,'/',comp_name));
    end

   for pair_num = 1:length(Ref_pair)
     pair_name = pai_name(pair_num+2).name;
    right_file = getdir(Right_pair{pair_num});
    fil_name = dir(Right_pair{pair_num});
    [reflagxx,refamp]=readsac(Ref_pair{pair_num});
         
    ref = refamp.*replace_lag';
    lagg = (1:length(refamp)) -mid-1;

% -----Ref stretching------ 
   str_coeff= linspace(-st_dura ,st_dura ,st_num+1);   

 for num =1:length(str_coeff)
    str_time = lagg./(1+str_coeff(num));
   % REF_str(num,:) = interp1(1:length(ref),ref,str_time+mid+1,'spline',0);  
   %   pre-filter is needed  map_coordinates is better
   REF_str(num,:) = double(ndimage.map_coordinates([ref zeros(1,length(ref))'] ...
        , [str_time+mid+1;zeros(1,length(ref))]));
 end
  %  cREF_str = foldlr(REF_str,crepl,crepr);
    sREF_str = allmove(REF_str,crepl,crepr);
% ----Ref Gauss4 fit Setup----
ft = fittype( 'gauss4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.Normalize = 'on';
opts.StartPoint = [0.028 -0.04 0.06 0.01 0.03 0.08 0.01 -0.2 0.1 0.008 0.2 0.1];

[XData, YData] = prepareCurveData( Rightlag, abs(refamp) );
[Fitresult, gof] = fit( XData, YData, ft, opts );

%  for Abnormal situations
Fitresult_names = coeffnames(Fitresult);
Confidence_interval = confint(Fitresult);
Abnormal_fit=abs(Confidence_interval(1,:)-Confidence_interval(2,:))>2;
for Abnormal_num = 1:3:length(Fitresult_names )
  if  abs(Fitresult.(Fitresult_names(Abnormal_num)))>0.1 ||...
  sum(Abnormal_fit(Abnormal_num:Abnormal_num+2))>3
  Fitresult.(Fitresult_names(Abnormal_num))=0;
  end
end

Symcenter = (Fitresult.b1*Fitresult.a1^2+Fitresult.b2*Fitresult.a2^2+ ...
 Fitresult.b3*Fitresult.a3^2+Fitresult.b4*Fitresult.a4^2)/(Fitresult.a1^2+...
 Fitresult.a2^2+Fitresult.a3^2+Fitresult.a4^2);

if isnan(Symcenter)
Symcenter = 0;
end
 

%% ----Prediction offset---
 for file_num = 1:length(right_file)
file_name{file_num} = fil_name(file_num+2).name;
      file_time=split(file_name(file_num),'.');
      time(file_num) = datenum(file_time(1));
     [~,rightamp] = readsac(right_file{file_num});  
     all_file(:,file_num) = rightamp ;

  %--- Right(days) Gauss4 fit Setup----
  
right_opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
right_opts.Display = 'Off';
right_opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
right_opts.Normalize = 'on';
right_opts.StartPoint = [max(rightamp) fix(mean(find(rightamp==max(rightamp))))/...
    sample_rate-sidewidth  std(rightamp)  0.025 0 std(rightamp) 0.01 -0.2 0.1 0.008 0.2 0.1];

   [xData, yData] = prepareCurveData(Rightlag, abs(rightamp));
   [fitresult, gof] = fit( xData, yData, 'Gauss4', right_opts );

   % for abnormal situation
fitresult_names = coeffnames(fitresult);
confidence_interval = confint(fitresult);
abnormal_fit=abs(confidence_interval(1,:)-confidence_interval(2,:))>2;
for abnormal_num = 1:3:length(fitresult_names )
  if  abs(fitresult.(fitresult_names(abnormal_num)))>0.1 ||...
  sum(abnormal_fit(abnormal_num:abnormal_num+2))>=3
  fitresult.(fitresult_names(abnormal_num))=0;
  end
end

   symcenter(file_num) = (fitresult.b1*fitresult.a1^2+fitresult.b2*fitresult.a2^2+ ...
   fitresult.b3*fitresult.a3^2+fitresult.b4*fitresult.a4^2)/(fitresult.a1^2+...
   fitresult.a2^2+fitresult.a3^2+fitresult.a4^2);

   if isnan(symcenter(file_num))
    symcenter(file_num) = 0;
   end

 end

%---------Get pair_offset------
 coffset = fix((Symcenter-symcenter)/0.05);   
 coffset(abs(coffset)>maxoffset)=0;
 coffset =  eliminate(coffset,maxdiff,dura_lag);
 [unique_values, ~, idx] = unique(coffset);
counts = histcounts(idx, length(unique_values));
 [~,index_count]=sort(counts,'descend');
if length(index_count) ==1 
     pair_offset = max(abs(unique_values(index_count)));
elseif length(index_count) == 2 && max(abs(unique_values(index_count(1:2))))<maxoffset
     pair_offset = max(abs(unique_values(index_count(1:2))));
else
  pair_offset = max(abs(unique_values(index_count(1:3))));
end


%% -----Prediction dvv-----

 for file_num = 1:length(right_file)

         rightamp = all_file(:,file_num); 
  right = allmove(rightamp',crepl,crepr,coffset(file_num));

%  coeffs = corr(right',cREF_str');
  coeffs = corr(right',sREF_str');
 max_index = fix(mean(find(coeffs == max(coeffs))));

%% ----loop for corrected offset----

   corre_dura =( 1:(2*(abs(pair_offset)+1)*4+1))-(abs(pair_offset)+1)*4-1;
   for  coffset_loop = 1:(2*(abs(pair_offset)+1)*4+1)
     corre_offse=coffset_loop-(abs(pair_offset)+1)*4-1;
   off_right = allmove(rightamp',crepl,crepr,corre_offse); 

%  corre_coeff(coffset_loop)=corr(cREF_str(max_index,:)',off_right'); %2->off_right'        
   corre_coeff(coffset_loop)=corr(sREF_str(max_index,:)',off_right'); %2->off_right'        
   
   end

   corre_index = find(corre_coeff==max(corre_coeff));  
   ccoffset(file_num) = fix(mean(corre_dura(corre_index)));
 end



%% ------ Correct DVV---

ccoffset(abs(ccoffset)>maxoffset)=0;
ccoffset =  eliminate(ccoffset,maxdiff,dura_lag);


 for file_num = 1:length(right_file)

        rightamp = all_file(:,file_num); 
   cright =  allmove(rightamp',crepl,crepr,ccoffset(file_num));

  ccoeffs = corr(cright',sREF_str');
%  ccoeffs = corr(cright',cREF_str');
  cmax_index = find(ccoeffs == max(ccoeffs)); 
  dvv(file_num) = -str_coeff(cmax_index);
  coeff(file_num) = max(ccoeffs);
 
 end    
 clear corre_coeff symcenter;

data_offset(offset_num,:) = ccoffset;
offset_num = offset_num+1;


writematrix([datenum(time)',dvv',coeff'], ...
       sprintf('%s/%s/%s.csv',patho,comp_name,pair_name));

  clear REF_str;
   end  
end

save(sprintf('%s/data_offset.mat',path_offset),'data_offset');
warning('on');

%% Function

function [paths] = getdir(path)
%'Test' -->'Test/MSN_test'
file=dir(path);
paths={};

 for i=3:length(file)
    right_path=strcat(path,'/',file(i).name);
    paths{i-2}=right_path;
 end
end


function [out] = allmove(one,crepl,crepr,offset)
% remain doble sides for chosen dura
  if nargin < 4 || offset==0     
    elseif   offset>0
    crepl = [zeros(1,offset) crepl(1:end-offset)];
    crepr = [zeros(1,offset) crepr(1:end-offset)];
    else
    crepl = [crepl(1-offset:end) zeros(1,-offset)];
    crepr = [crepr(1-offset:end) zeros(1,-offset)];
  end  
right = ismember(crepr+crepl,1);
out = one(:,right);

end


function [out] = foldlr(one,crepl,crepr,offset)
%-----move Fold average-----
    if nargin < 4 || offset==0     
    elseif   offset>0
    crepl = [zeros(1,offset) crepl(1:end-offset)];
    crepr = [zeros(1,offset) crepr(1:end-offset)];
    else
    crepl = [crepl(1-offset:end) zeros(1,-offset)];
    crepr = [crepr(1-offset:end) zeros(1,-offset)];
    end  

   out = (fliplr(one(:,crepl>0))+one(:,crepr>0))/2; 

end



function    [offset] =  eliminate(offset,maxdiff,dura_lag)

 %  cut oscillating offset (fake clockerror)

 % obviously, there are some abnormal situations 
 % I just ignore some of these  because of if if if if  
 %   Gauss maybe more easier  --> eliminate the edge
 
 % [coffset;oscillate 0;bool 0 0];

  diff_offset = abs(diff(offset));  
  oscillate = diff_offset>maxdiff;
   bool = diff(oscillate);osc_node=1;   %fir_tran=0;

   start_=find(bool==1); end_=find(bool==-1);
if isempty(start_) && isempty(end_)
elseif  length(start_)==length(end_) && start_(1)<end_(1) 

       for   i = 1:length(start_) 
           if   any(end_(i)-start_(i) == 1)
               offset(start_(i)+2) = 0;     
           else
               offset(start_(i)+2:end_(i)) = 0;     
           end
       end

elseif   start_==length(offset)-2 
      offset(end) =0;
elseif     end_ == 1
      offset(1) =0;

else
   for bool_num =1 :length(bool)-1
     if bool(bool_num)  
       if bool(bool_num)>0 &&(osc_node == 1|| ...
             bool_num>dura_lag && ~any(bool(bool_num-dura_lag:bool_num-1)) && ~bool(bool_num+1))
          trans(osc_node)=bool_num+1;
         osc_node=osc_node+1;
       elseif  bool_num>bool_num && ~any(bool(bool_num-dura_lag:bool_num-1)) && ~bool(bool_num+1)...
           ||  length(bool)>=bool_num+dura_lag && ~any(bool(bool_num+1:bool_num+dura_lag))...

          if osc_node == 1 
             aaa=find(bool<0);     
             trans(1)= aaa(1);       trans(2)=bool_num+1; 
          end

         trans(osc_node)=bool_num+1;
         osc_node=osc_node+1;
       elseif length(bool)<bool_num+dura_lag
          trans(osc_node)=bool_num+1; break;
       end
     end
   end

 
   if length(trans)==1
  if length(offset)/2>trans
     offset(1:trans)=0;
  else
       offset(trans:end)=0;
  end
 elseif ~oscillate(1)
    for zer=1:2:length(trans)
      if zer ==length(trans) && zer~=1
    offset(trans(zer):end) = 0;break;
      end
     offset(trans(zer)+1:trans(zer+1)-1) = 0;
    end
  else
    for zer=1:2:length(trans)
      if zer ==1 
    offset(1:trans(zer)) = 0; 
         if zer ==length(trans)-1
    offset(trans(zer+1):end) = 0;break;
         end
       continue;
      end                                                                                   

    offset(trans(zer-1)+1:trans(zer)-1) = 0;

     if zer ==length(trans)-1
    offset(trans(zer+1):end) = 0;break;
     end

    end
  end

end

end
