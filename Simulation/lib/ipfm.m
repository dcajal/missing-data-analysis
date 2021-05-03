function [ mt , d_HR ] = ipfm ( tm, ids, t , ord )
% Estimation of instantaneous heart rate based on the IPFM model
% sp        =	IPFM (tn,ids,ord)
% tn        =	normal beat ocurrence time series (in seconds)
% ids       =	beat indices following each incidence
% t         =   vector of time with query points for interpolation
% ord       =	order of the spline for the interpolation (default, 14)
% d_HR(t)	=	(1+m(t))/T interpolated with splines of order 'ord'

% (c) <Javier Mateo>, 7-Dic-2001. 
% Modifications <Eduardo Gil>,<Raquel Bailón>, 30-Dic-2009
% 'Simplified' and 'documented' version. Previous name: mti  <Pablo Armañac> 2019
% 
%

% s represents the phase shift in the integration process of the IPFM
if nargin<4
    ord =   5  ;
end

if isempty(tm)
    error('The Normal Beat Ocurrence time series is empty');
end

% % % % % % % % % % % % % % % % % % % % % 
tic
% % % % % % % % % % % % % % % % % % % % % 

tm  =   tm(:)   ;
ids =   ids(:)  ;
s   =   []      ;

MG  =   2       ;
k   =   10      ;

dt1 =   median( diff(tm(1:9)));
dt2 =   median( diff(tm(end-8:end)));
tt  =   [ tm(1)-(k:-1:1)'*dt1 ; tm ; tm(end)+(1:k)'*dt2 ] - tm(1);

if ~isempty(ids)
    try
        [s,tex,x]   =   calcjump ( tt , ids+k , 4 , MG )    ;
    catch
        disp('Timeout reached. mt set as empty.');
        mt = [];
        d_HR = [];
        return
    end
else
    s           =   []                                  ; 
    tex         =   tt                                  ; 
    x           =   (1:length(tt))'                     ;
end

d_HR	=   fnder(  spapi( ord , tex+tm(1) , x )  );

mt      =	spval(d_HR,t);

% Zeros not allowed
if mt(1) == 0
    mt(1:find(mt>0,1)) = mt(find(mt>0,1));
end
if mt(end) == 0
    mt(find(mt>0,1,'last'):end) = mt(find(mt>0,1,'last'));
end

end


function [ s , tex , x ] = calcjump ( tt , ids , ord , MG )
% method for determining these a-priori unknown s and T quantities required to estimate ht(tk)

% the beat timing after the ectopic beat is not a delayed version of the beat timing that would have occurred if the ectopic beat had not occurred. 
% Dealing with ectopic beats correction by introducing a constant delay in the beat timing cannot then correctly recover the beat sequence. 
% At best, the delay could only absorb the spike at the ectopic beat position but this will result in a phase shift in the HT or other HRV signals 
% before and after of the ectopic beat, and consequently, the spectra will continue to be corrupted.

hp          =   diff(tt)                            ;
[xn,yn]     =   fillgap(tt,hp,ids,MG)               ;
sp          =   spapi(ord,xn,yn)                    ;
tnew        =   calctnew(sp,tt(ids-1),tt(ids),ids)  ;

tex{1}      =   tt(1:ids(1)-1)                      ;
x{1}        =   (1:ids(1)-1)'                       ;
spx(1)      =   spapi(ord,tex{1},x{1})              ;
ids         =   [ids;length(tt)+1]                  ;

for i=1:length(ids)-1
    tex{i+1}    =   [tnew{i}';tt(ids(i):ids(i+1)-1)]        ;
    x{i+1}      =   (ids(i)-length(tnew{i}):ids(i+1)-1)'    ;
    spx(i+1)    =   spapi(ord,tex{i+1},x{i+1})              ;
    s(i)        =   spval(spx(i),tex{i+1}(1))-x{i+1}(1)     ;
    x{i+1}      =   x{i+1}+sum(s)                           ;
end

tex         =   cat(1,tex{:});
x           =   cat(1,x{:});

% A safe interpretation of the jump 's' is as follows:
% s = 1 means that something occurred (probably a false negative or a ventricular ectopic beat) that did not reset the process of beats generation  [see Fig. 1(b), (d)]
% s < 1 means that something occurred (probably a supraventricular ectopic beat) that reset the process of beats generation prematurely             [see Fig. 1(c)]
% s > 1 means that several consecutive normal beats were missed due to ectopic beats or false negatives                                             [see Fig. 1(e), (f)]
% False positives are simply corrected by deleting them                                                                                             [see Fig. 1(a)]

end

function [ xn , yn ] = fillgap ( xx , yy , ids , MG )
%length(xx)=length(yy)+1
%Ej: xx=tn, yy=diff(tn) o yy=60./diff(tn)
%Ojo: xn e yn estan sin ordenar

idv     =   setdiff((1:length(yy))',ids-1)  ;
xn      =   xx(idv+1)                       ;
yn      =   yy(idv)                         ;
dx      =   diff(xn)                        ;
mfdx    =   MG*medfilt1(dx,21)              ;
idg     =   find(dx>mfdx)                   ;

for i=1:length(idg)
    
    id      =   idg(i)                      ;
    if id >= length(mfdx)
        delta = mfdx(end);
    else
        delta   =   min(mfdx(id:id+1))          ;
    end
    
    xi      =   xn(id)                      ;
    xf      =   xn(id+1)                    ;
    
    yi      =   yn(id)                      ;
    yf      =   yn(id+1)                    ;
    
    stp     =   ceil( (xf-xi)/delta )       ;
    cc      =   (1:stp-1)'./stp             ;
    
    xn      =   [xn; xi+cc*(xf-xi) ]         ;
    yn      =   [yn; yi+cc*(yf-yi) ]         ;
    
end

end

function tnew = calctnew ( sp , ta , tb , ids )

% Instantaneous Heart Rate and HT Determination: 
% The starting point of this problem is the known occurrence times of the normal beats before, tk = t(k) with k<ke, and after the ectopic beat, tk = t(k-1+s) with k>ke. 
% Using these values, the samples of the HP signal, hp(t{k}) = t{k}-t{k-1}, can be calculated, taking into account that all values obtained will be correct except for those involving the ectopic beat. 
% So, it can be seen, for the beats before the ectopic beat, the values of the HP signal are hp(t(k)) = t(k)-t(k-1) with k<ke, 
% and for the beats after the ectopic beat hp(t(k-1+s)) = t(k-1+s)-t(k+s-2) with k>ke+1.
% This incorrect value must be rejected before any further assessment. 
% Using only the correct values of allows a continuous time estimation by means of spline interpolation.

% Once estimated hp_est(t), the virtual beat times t(ke-n+s), being t(ke-n+s) the hypothetical nth normal beat positions backward extended prior to the first normal beat (t(ke+s)) 
% after the ectopic beat, can be calculated as t(ke-n+s) = t(ke-n+1+s) - hp_est( t(ke-n+1+s) )

for i=1:length(ids)
    j           =   1                               ;
    tnew{i}(1)  =   soltb(sp,tb(i))                 ;
    d           =   ta(i)-tnew{i}(1)                ;
    while d<0
        tnew{i}(j+1)    =   soltb(sp,tnew{i}(j))    ;
        j               =   j+1                     ;
        d               =   ta(i)-tnew{i}(j)        ;
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        if toc > 1*10
            disp('% % % % % % % % % % % % % % % % % % % % % % % %');
            error('m(t) Estimation FAILED Succesfully');
            disp('% % % % % % % % % % % % % % % % % % % % % % % %');
            break
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
    end
    tnew{i}     =   fliplr( tnew{i} )               ;
end

end


function tb = soltb ( hpsp , tk )

tb  =   tk - spval( hpsp , tk ) ;

end



% function tnew=calctnew(sp,ta,tb,ids)
% % With an attempt to break when the code stucks in an infinite loop
% for i=1:length(ids)
%     j=1;
%     tnew{i}(1)=soltb(sp,tb(i));
%     d=ta(i)-tnew{i}(1);
%     
%     time0 = tic;
%     timeLimit = 15; % 1 hour == 3600 seconds
%     
%     while d<0
%         tnew{i}(j+1)=soltb(sp,tnew{i}(j));
%         j=j+1;
%         d=ta(i)-tnew{i}(j);
%         if toc(time0)>timeLimit
%             disp('se cuelga en el mti');
%             break
%         end
%     end
%     
%     tnew{i}=fliplr(tnew{i});
%     
% end
% end
