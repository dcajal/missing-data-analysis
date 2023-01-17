function [sp] = mti(tn,ids,ord)
% Estimation of instantaneous heart rate based on the IPFM model
% sp  = mti(tn,ids,ord)
% tn  = normal beat ocurrence time series (in seconds)
% ids = beat indices following each incidence
% ord = order of the spline for the interpolation (default, 14)
% sp  = (1+m(t))/T interpolated with splines of order ord 
% (c) Javier Mateo, 7-Dic-2001. Modified by Eduardo Gil and Raquel Bailón,
% 30-Dic-2009
%
if nargin<3, ord=14; end
MG=2;
if isempty(tn)
    sp=[];s=[];
    return
end
tn=tn(:);ids=ids(:);s=[];

k=10;
dt1=median(diff(tn(1:9)));
dt2=median(diff(tn(end-8:end)));
tt=[tn(1)-(k:-1:1)'*dt1;tn;tn(end)+(1:k)'*dt2]-tn(1);

if ~isempty(ids)
    [s,tex,x] = calcjump(tt,ids+k,4,MG); 
else
    s=[]; tex=tt; x=(1:length(tt))';
end

sp=fnder(spapi(ord,tex+tn(1),x));


% %%
% labels{1} = 't_mti';
% desc = labels;
% units = {[]};
% for i = 1:size(labels,1)
%     units(i,1) = {'ms'};
% end
% t_mti = tex+tn(1);
% t_mti(1:10) = [];
% t_mti(end-10:end) = [];
% ID_nums = 1:length(t_mti);
%  eval(['v.' labels{i} '.pos = t_mti*1000;']);
%  eval(['v.' labels{i} '.info{1} = zeros(length(v.' labels{i} '.pos),3);']);
%  eval(['v.' labels{i} '.info{1}(:,1) = ID_nums;']);
%  eval(['v.' labels{i} '.info{2}(1:length(v.' labels{i} '.pos)) = {''''};']);
% %%
% global prmSignal
% writeannot_mat(prmSignal.annotdir, prmSignal.ecgannot(1:end-4),labels',desc',units',v); %#ok<NODEF>

function [s,tex,x] = calcjump(tt,ids,ord,MG)
hp=diff(tt);
[xn,yn]=fillgap(tt,hp,ids,MG);
sp=spapi(ord,xn,yn);
tnew=calctnew(sp,tt(ids-1),tt(ids),ids);
tex{1}=tt(1:ids(1)-1);
x{1}=(1:ids(1)-1)';
spx(1)=spapi(ord,tex{1},x{1});
ids=[ids;length(tt)+1];
for i=1:length(ids)-1
    tex{i+1}=[tnew{i}';tt(ids(i):ids(i+1)-1)];
    x{i+1}=(ids(i)-length(tnew{i}):ids(i+1)-1)';
    spx(i+1)=spapi(ord,tex{i+1},x{i+1});
    s(i)=spval(spx(i),tex{i+1}(1))-x{i+1}(1);
    x{i+1}=x{i+1}+sum(s); 
end
tex=cat(1,tex{:});
x=cat(1,x{:});

function [xn,yn]=fillgap(xx,yy,ids,MG)
%length(xx)=length(yy)+1
%Ej: xx=tn, yy=diff(tn) o yy=60./diff(tn)
%Ojo: xn e yn estan sin ordenar
%
idv=setdiff([1:length(yy)]',ids-1);
xn=xx(idv+1);yn=yy(idv);
dx=diff(xn);
mfdx=MG*medfilt1(dx,21);
idg=find(dx>mfdx);
for i=1:length(idg)
    id=idg(i);
    if id >= length(mfdx)
        delta = mfdx(end);
    else
        delta = min(mfdx(id:id+1));
    end
    xi=xn(id);xf=xn(id+1);
    yi=yn(id);yf=yn(id+1);
    stp=ceil((xf-xi)/delta);
    cc=(1:stp-1)'./stp;
    xn=[xn;xi+cc*(xf-xi)];
    yn=[yn;yi+cc*(yf-yi)];
end

function tnew=calctnew(sp,ta,tb,ids);
for i=1:length(ids)
    j=1;
    tnew{i}(1)=soltb(sp,tb(i));
    d=ta(i)-tnew{i}(1);
    while d<0
        tnew{i}(j+1)=soltb(sp,tnew{i}(j));
        j=j+1;
        d=ta(i)-tnew{i}(j);
    end
    tnew{i}=fliplr(tnew{i});
end

function ta=solta(hpsp,tk)
ta=tk;
err=1;
while err>1e-6
   tanew=spval(hpsp,ta)+tk;
   err=abs(ta-tanew);
   ta=tanew;
end

function tb=soltb(hpsp,tk)
tb=tk-spval(hpsp,tk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp = spapi(knots,x,y)
%SPAPI Spline interpolation.
%
%   SPAPI(KNOTS,X,Y)      (with both KNOTS and X vectors)
%   returns the spline  f  (if any) of order
%         k := length(KNOTS) - length(X)
%   with knot sequence KNOTS for which  Y(:,j) = f(X(j)) , all j.
%   This is taken in the osculatory sense in case some X are repeated,
%   i.e., in the sense that  D^m(i) f(X(i)) = Y(:,i)  in case the X are
%   in nondecreasing order, with  m(i) := #{ j<i : X(j) = X(i) }.
%   For this, it is important to realize the following: if the input X is
%   not nondecreasing, it will be reordered to make it so, and the input Y 
%   will be reordered in the same way. See the cubic Hermite interpolation 
%   example below.
%
%   SPAPI(K,X,Y)          (with K a positive integer)
%   uses APTKNT to obtain from K and X an apt knot sequence. Precisely,
%   it has the same effect as  SPAPI(APTKNT(sort(X),K),X,Y).
%
%   SPAPI({KNOTS1,...,KNOTSm},{X1,...,Xm},Y)
%   returns in SP the m-variate tensor-product spline of coordinate order
%      ki := length(KNOTSi) - length(Xi)
%   with knot sequence KNOTSi in the i-th variable, i=1,...,m, for which
%   Y(:,j1,...,jm) = f(X1(j1),...,Xm(jm)),  all j := (j1,...,jm) 
%   As in the univariate case, KNOTSi may also be a positive integer, in which
%   case the i-th knot sequence is obtained from Xi via APTKNT.
%   Note the possibility of interpolating to vector-valued data. However,
%   in contrast to the univariate case, if the data to be interpolated are
%   scalar-valued, then the input array Y is permitted to be m-dimensional,
%   in which case
%   Y(j1,...,jm) = f(X1(j1),...,Xm(jm)),  all j := (j1,...,jm) 
%   Multiplicities in the sequences Xi lead to osculatory interpolation just
%   as in the univariate case.
%
%   For example, if the points in the vector  t  are all distinct, then
%
%      sp = spapi(augknt(t,4,2),[t t],[cos(t) -sin(t)]);
%
%   provides the C^1 piecewise cubic Hermite interpolant to the function 
%   f(x) = cos(x)  at the points in  t . If matching of slopes is only
%   required at some subsequence  s  of  t  but that includes the leftmost
%   and the rightmost point in  t , one would use instead
%
%      sp = spapi( augknt([t s],4), [t s], [cos(t) -sin(s)] );
%
%   or even just
%
%      sp = spapi( 4, [t s], [cos(t) -sin(s)] );
%
%   and the last works even if s fails to include the extreme points of t,
%   and produces even a C^2 piecewise cubic interpolant to these Hermite data.
%
%   As another example,
%
%      sp = spapi( {[0 0 1 1],[0 0 1 1]}, {[0 1],[0 1]}, [0 0;0 1] );
%
%   constructs the bilinear interpolant to values at the corners of a square,
%   as would the statement
%
%      sp = spapi({2,2},{[0 1],[0 1]},[0 0;0 1]);
%
%   As a final example, the statements
%
%      x = -2:.5:2; y=-1:.25:1; [xx, yy] = ndgrid(x,y); 
%      z = exp(-(xx.^2+yy.^2)); 
%      sp = spapi({3,4},{x,y},z);
%      fnplt(sp)
%
%   produce the picture of an interpolant (piecewise quadratic in x, 
%   piecewise cubic in y) to a bivariate function.
%   Use of MESHGRID instead of NDGRID here would produce an error.
%
%   See also SPAPS, SPAP2.

%   Carl de Boor 24 dec 89
%   cb :  22 apr 95
%   cb :  9 may '95 (use .' instead of ')
%   cb :  9 mar '96 (flip  X  in case it has more than one row)
%   mathworks: 19mar96 (clean up error(...))
%   cb : 12sep96 (string argument for SPCOL)
%   cb : 18oct97 (include tensor product)
%   cb : 28nov98 (handle sparse Y, add bivariate example, mention NDGRID)
%   cb : 14feb99 (with KNOTS a positive integer, take it as the order and use
%                 APTKNT to construct a knot sequence apt for the given X)
%   cb : 21nov99, 20feb00 (improve the help, permit k=1)
%   cb : 10apr00 (introduce optional third argument for multidimensional case;
%                 also handle just one data point in case k=1;
%                 also use the actual order used in APTKNT in case the user
%                 specified the order)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.15 $

if iscell(knots) % gridded data are to be interpolated by tensor product splines
   if ~iscell(x)
      error('If KNOTS is a cell-array, then also X must be one.'), end
   m = length(knots);
   if m~=length(x)
      error('If KNOTS is a cell-array, then X must be one of the same length.')
   end
   sizey = size(y);
   switch length(sizey)
     case m  % grid values of a scalar-valued function
        if issparse(y), y = full(y); end 
        sizey = [1 sizey]; y = reshape(y, sizey); 
     case m+1
     otherwise
        error(['If KNOTS is a cell-array of length m, then Y must be an ', ...
               'm- or (m+1)-dimensional array.'])
   end
   
   v = y; sizev = sizey;
   for i=m:-1:1   % carry out coordinatewise interpolation
      temp = spapi1(knots{i},x{i}, reshape(v,prod(sizev(1:m)),sizev(m+1)));
      v = reshape(spbrk(temp,'c'), sizev);
      if length(knots{i})==1, knots{i} = spbrk(temp,'knots'); end
      if m>1
         v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
      end
   end
   % At this point, V contains the tensor-product B-spline coefficients;
   % also, the various knot sequences will have been checked. 
   % It remains to return information:
   sp = spmak(knots, v, sizev);

else             % univariate spline interpolation
   sp = spapi1(knots,x,y);
end

function sp = spapi1(knots,x,y)
%SPAPI1 univariate spline interpolation

if length(knots)==1 % the order is being specified; get the knots via APTKNT
   k = knots;
   if k~=fix(k), error('The order must be a natural number.'), end
   if k<1, error('The order must be at least 1.'), end
   if k==1&length(x)==1
      knots = [x, x+1];
   else
      knots = x; if any(diff(knots)<0), knots = sort(knots); end
      [knots,k] = aptknt(knots,k);
   end
else

   if ~isempty(find(diff(knots)<0))
      error('The knot sequence should be nondecreasing.')
   end
   
   n = length(x); npk = length(knots); k = npk-n;
   if k<1
     error(sprintf(['The number of data points, %.0f, should be less than\n',...
             'the number, %.0f, of knots.'],n, npk));
   end
end

% It is assumed that  y  has the same sense as  x  which is important in case
% the entries of  y  are vectors rather than scalars
[rx,cx] = size(x); [ry,cy] = size(y);
if rx>1
   if ry~=rx
     error(sprintf(['Since  x  has more than one row,\n',...
            '   y  should have the same number of rows as  x .']));
   else, x = x.'; y = y.'; end
else,
   if cy~=cx error('Number of data sites and data values should match.')
   end
end

%  sort the given abscissae.
[x,index] = sort(x); y = y(:,index);

%  Generate the collocation matrix and divide it into the possibly reordered
%  sequence of given ordinates to generate the B-spline coefficients of the
%  interpolant, then put it all together into SP.

sp = spmak(knots,slvblk(spcol(knots,k,x,'slvblk'),y.').');

function colloc = spcol(knots,k,tau,arg1,arg2)
%SPCOL B-spline collocation matrix.
%
%   SPCOL(KNOTS,K,TAU)  returns the matrix   (D^m(i)B_j(TAU_i)) , with B_j
%   the j-th B-spline of order K for the knot sequence KNOTS, TAU a
%   sequence of points, both assumed to be  n o n d e c r e a s i n g ,
%   and   m(i) := #{ j<i : TAU(j) = TAU(i) }.
%
%   If one of the optional arguments is a string with the same first two
%   characters as 'slvblk', the matrix is returned in the almost block-
%   diagonal format (specialized for splines) required  by SLVBLK.m .
%   (For backward compatibility, this also happens if the fourth argument is
%   present and is not a string).
%
%   If one of the optional arguments is a string with the same first two
%   characters as 'sparse', the matrix is returned as a sparse matrix.
%
%   If one of the optional arguments is a string with the same first character
%   as 'noderiv', then multiplicities are ignored, i.e., m(i) := 1 for all i.
%   (For backward compatibility, this also happens if the fifth argument is
%   present and is not a string).
%
%   The recurrence relations are used to generate the entries of the
%   matrix.
%
%   For example, the statement
%
%      x = linspace(t(1),t(end),101); c = spcol(t,3,x);
%
%   provides, in c(:,j), a fine sequence of values of the j-th quadratic 
%   B-spline for the given knot sequence t .
%
%   See also SLVBLK, SPARSE, SPAPI, SPAP2, BSPLINE.

%   Carl de Boor 13 may 89
%   cb : 19 apr 1993; correct output for case  n < 1
%   cb :  7 jun 1994; add option to ignore multiplicities
%   cb :  9 may '95 (use .' instead of ')
%   cb :  8 dec '95 (update help, add check on TAU and on M)
%   cb : 22 may '96 (make optional arguments more flexible)
%   cb : 11sep96 (add the option of returning a sparse matrix)
%   cb : 5 oct 97 (replace use of ANS, to help compilation)
%   cb : 19may98 (standardize the help)
%   cb : 28feb99 (use repmat; accept TAU as row or column)
%   cb : 24aug99 (correct mispelling in help)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.12 $

if ~isempty(find(diff(knots)<0))
   error('The knot sequence KNOTS should be nondecreasing.')
end
if ~isempty(find(diff(tau)<0))
   error('The point sequence TAU should be nondecreasing.')
end

%  Compute the number  n  of B-splines of order K supported by the given
%  knot sequence and return empty matrix in case there aren't any.

npk=length(knots); n=npk-k;
if n<1 fprintf('There are no B-splines for the given input.\n')
   colloc = []; return
end

% Settle the options:
slvblk=0; noderiv=0;
if nargin>3
   if isstr(arg1)
      if     arg1(1)=='s', slvblk =1;
         if length(arg1)>1, if arg1(2)=='p', slvblk = 2; end, end
      elseif arg1(1)=='n', noderiv=1;
      else error(['The second optional argument should be ''sl'', ''sp'', ',...
      '''n'', or a number.'])
      end
   else slvblk=1;
   end
end
if nargin>4
   if isstr(arg2)
      if     arg2(1)=='s', slvblk =1;
         if length(arg2)>2, if arg2(2)=='p', slvblk = 2; end, end
      elseif arg2(1)=='n', noderiv=1;
      else error(['The second optional argument should be ''sl'', ''sp'', ',...
      '''n'', or a number.'])
      end
   else noderiv=1;
   end
end

% If  NODERIV==0, remove all multiplicities from TAU and generate repetitions
% of rows instead of rows containing values of successive derivatives.
nrows = length(tau); tau = reshape(tau,1,nrows);
if noderiv
   index = 1:nrows; m = ones(1,nrows); nd = 1; pts = tau;
else
   index = [1 find(diff(tau)>0)+1];
   m = diff([index nrows+1]); nd = max(m);
   if nd>k
      error(sprintf('Point multiplicity should not exceed the given order %g.',k));
   end
   pts = tau(index);
end

%set some abbreviations
km1 = k-1;

%  augment knot sequence to provide a K-fold knot at each end, in order to avoid
% struggles near ends of basic interval,  [KNOTS(1) .. KNOTS(npk)] .
% The resulting additional B-splines, if any, will NOT appear in the output.

[augknot,addl] = augknt(knots,k); naug = length(augknot)-k;
pts = pts(:); augknot = augknot(:);

%  For each  i , determine  savl(i)  so that  K <= savl(i) < naug+1 , and,
% within that restriction,
%        augknot(savl(i)) <= pts(i) < augknot(savl(i)+1) .

savl = max([sorted(augknot(1:naug),pts); k*ones(1,length(pts))]);

b = zeros(nrows,k);

% first do those without derivatives
index1 = find(m==1);
if ~isempty(index1)
   pt1s = pts(index1); savls = savl(index1); lpt1 = length(index1);
   % initialize the  b  array.
   lpt1s = index(index1); b(lpt1s,1) = ones(lpt1,1);

   % run the recurrence simultaneously for all  pt1(i) .
   for j=1:km1
      saved = zeros(lpt1,1);
      for r=1:j
         tr = augknot(savls+r)-pt1s;
         tl = pt1s-augknot(savls+r-j);
         term = b(lpt1s,r)./(tr+tl);
         b(lpt1s,r) = saved+tr.*term;
         saved = tl.*term;
      end
      b(lpt1s,j+1) = saved;
   end
end

% then do those with derivatives, if any:
if nd>1
   indexm=find(m>1);ptss=pts(indexm);savls=savl(indexm);lpts=length(indexm);
   % initialize the  bb  array.
   temp = [1 zeros(1,km1)]; bb = temp(ones(nd*lpts,1),:);
   lptss = nd*[1:lpts];

   % run the recurrence simultaneously for all  pts(i) .
   % First, bring it up to the intended level:
   for j=1:k-nd
      saved = zeros(lpts,1);
      for r=1:j
         tr = augknot(savls+r)-ptss;
         tl = ptss-augknot(savls+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptss,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptss,j+1) = saved;
   end

   % save the B-spline values in successive blocks in  bb .

   for jj=1:nd-1
      j = k-nd+jj; saved = zeros(lpts,1); lptsn = lptss-1;
      for r=1:j
         tr = augknot(savls+r)-ptss;
         tl = ptss-augknot(savls+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptsn,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptsn,j+1) = saved; lptss = lptsn;
   end

   % now use the fact that derivative values can be obtained by differencing:

   for jj=nd-1:-1:1
      j = k-jj;
      temp = [jj:nd-1].'*ones(1,lpts)+ones(nd-jj,1)*lptsn; lptss=temp(:);
      for r=j:-1:1
         temp = ones(nd-jj,1)*(augknot(savls+r)-augknot(savls+r-j)).'/j;
         bb(lptss,r) = -bb(lptss,r)./temp(:);
         bb(lptss,r+1) = bb(lptss,r+1) - bb(lptss,r);
      end
   end

   % finally, combine appropriately with  b  by interspersing the multiple
   % point conditions appropriately:
   dtau = diff([tau(1)-1 tau(:).' tau(nrows)+1]);
   index=find(min(dtau(2:nrows+1),dtau(1:nrows))==0); % Determines all rows
                                                    % involving multiple tau.
   dtau=diff(tau(index));index2=find(dtau>0)+1;     % We need to make sure to
   index3=[1 (dtau==0)];                            % skip unwanted derivs:
   if ~isempty(index2)
             index3(index2)=1+nd-m(indexm(1:length(indexm)-1));end
   b(index,:)=bb(cumsum(index3),:);

   % ... and appropriately enlarge  savl
   index = cumsum([1 (diff(tau)>0)]);
   savl = savl(index);
end

% Finally, zero out all rows of  b  corresponding to TAU outside the basic
% interval,  [knots(1) .. knots(npk)] .

index = find(tau<knots(1)|tau>knots(npk));
if ~isempty(index)
   temp = repmat([1-nd:0].',1,length(index))+repmat(nd*index(:).',nd,1);
   b(temp(:),:) = zeros(nd*length(index),k);
end

% The first B-spline of interest begins at KNOTS(1), i.e., at  augknot(1+addl)
% (since  augknot's  first knot has exact multiplicity K). If  addl<0 ,
% this refers to a nonexistent index and means that the first  -addl  columns
% of the collocation matrix should be trivial.  This we manage by setting
savl = savl+max(0,-addl);

if slvblk     % return the collocation matrix in almost block diagonal form.
              % For this, make the blocks out of the entries with the same
              %  SAVL(i) , with  LAST  computed from the differences.
   % There are two issues, the change in the B-splines considered because of
   % the use of  AUGKNOT  instead of  KNOTS , and the possible drop of B-splines
   % because the extreme  TAU  fail to involve the extreme knot intervals.

   % SAVL(j) is the index in  AUGKNOT  of the left knot for  TAU(j) , hence the
   % corresponding row involves  B-splines to index  savl(j) wrto augknot, i.e.,
   % B-splines to index  savl(j)-addl  wrto  KNOTS.
   % Those with negative index are removed by cutting out their columns (i.e.,
   % shifting the blocks in which they lie appropriately). Those with index
   % greater than  n  will be ignored because of  last .

   if addl>0   % if B-splines were added on the left, remove them now:
      width = km1+k;cc = zeros(nrows*width,1);
      index = min(k,savl-addl); 
      temp = +repmat(nrows*[0:km1],nrows,1);
    cc(repmat(([1-nrows:0]+nrows*index).',1,k)+repmat(nrows*[0:km1],nrows,1))=b;
      b(:)=cc(repmat([1-nrows:0].',1,k)+repmat(nrows*(k+[0:km1]),nrows,1));
      savl=savl+k-index;
   end
   ds=(diff(savl));
   index=[0 find(ds>0) nrows];
   rows=diff(index);
   nb=length(index)-1;
   last=ds(index(2:nb));
   if addl<0  nb=nb+1; rows=[0 rows]; last=[-addl last]; end
   addr=naug-n-addl;
   if addr<0  nb=nb+1; rows=[rows 0]; last=[last -addr]; end
   if slvblk==1
      colloc=[41 nb rows k last n-sum(last) b(:).'];
   else   % return the equivalent sparse matrix (cf BKBRK)
      nr = (1:nrows).'; nc = 1:k; nrnc = nrows*k;
      ncc = zeros(1,nrows); ncc(1+cumsum(rows(1:(nb-1)))) = last;
      ncc = reshape(cumsum(ncc),nrows,1);
      ijs = [reshape(nr(:,ones(1,k)),nrnc,1), ...
           reshape( ncc(:,ones(1,k))+nc(ones(nrows,1),:) , nrnc,1), ...
           reshape(b,nrnc,1)];
      index = find(ijs(:,2)>n);
      if ~isempty(index), ijs(index,:) = []; end
      colloc = sparse(ijs(:,1),ijs(:,2),ijs(:,3),nrows,n);
   end
else          % return the collocation matrix in standard matrix form
   width = max([n,naug])+km1+km1;
   cc = zeros(nrows*width,1);
   cc(repmat([1-nrows:0].',1,k)+ ...
              repmat(nrows*savl.',1,k)+repmat(nrows*[-km1:0],nrows,1))=b;
   % (This uses the fact that, for a column vector  v  and a matrix  A ,
   %  v(A)(i,j)=v(A(i,j)), all i,j.)
   colloc = reshape(cc(repmat([1-nrows:0].',1,n) + ...
                    repmat(nrows*(max(0,addl)+[1:n]),nrows,1)), nrows,n);
end

function v = spval(sp,x,left)
%SPVAL Evaluate function in B-form.
%
%   SPVAL(SP,X)  returns the value at X of the function whose B-form is in SP.
%
%   SPVAL(PP,X,'l<anything>')  takes the function to be left-continuous.
%   If the function in SP is m-variate and the third argument is actually an
%   m-cell, LEFT say, then, continuity from the left is enforced in the i-th 
%   variable if  LEFT{i}(1) = 'l'.  
%
%   The output is a matrix of size [d*m,n] if the function in SP is
%   univariate and d-vector valued, and [m,n] == size(X) .
%
%   If the function in SP is m-variate with m>1 and d-vector valued, then
%
%                         [d,n],         if X is of size [m,n]
%   the output is of size [d,n1,...,nm], if d>1  and X is {X1,...,Xm} 
%                         [n1,...,nm],   if d==1 and X is {X1,...,Xm} 
%
%   See also FNVAL, PPUAL, RSVAL, PPVAL.

%   Carl de Boor 23 may 89
%   cb / 18 apr 93: fix case k = 1; prettify
%   cb / 05 nov 94: prettify
%   cb :  9 may '95 (use .' instead of ')
%   cb : 15sep96 (handle empty X, return matri of same size,
%                 edited help accordingly, simplify reordering,
%                 do replication without multiplication.)
%   cb : 5 oct 97 (replace use of ANS, to help compilation)
%   cb : 18oct97 (include tensor product)
%   cb : 10may98 (avoid replication by multiplication; standardize help)
%   cb : 11dec98 (use repmat, improve warning)
%   cb : 04may00 (introduce optional third argument for getting left continuity)
%   cb : 23jun00 (enforce claim that only a LEFT of the form 'l...' will cause
%                 left continuity; streamline determination of correct knot
%                 interval)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.11 $


if ~isstruct(sp)
   error('SP must be a structure.'), end
if iscell(sp.knots)  % we are dealing with a tensor product spline

   [t,a,n,k,d] = spbrk(sp); m = length(t);

   if nargin>2 % set up left appropriately
      if ~iscell(left)
         temp = left; left = cell(1,m); [left{:}] = deal(temp);
      end
   else
      left = cell(1,m);
   end

   if iscell(x)  % evaluation on a mesh

      v = a; sizev = size(v); nsizev = zeros(1,m);
      if length(sizev)~=m+1, error('Information in SP is inconsistent.'), end

      for i=m:-1:1
         nsizev(i) = length(x{i}(:));
         v = reshape(...
         spval1(spmak(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x{i},left{i}),   [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         v = reshape(v,[d,nsizev]);
      else
         v = reshape(v,nsizev);
      end
 
   else          % evaluation at scattered points;
                 % this will eventually be done directly here.
  fprintf('Warning (from SPVAL): Converted SP to ppform and then used PPUAL.\n')
      v = ppual(sp2pp(sp),x);
   end

else                 % we are dealing with a univariate spline
   if nargin<3, left = []; end
   v = spval1(sp,x,left);
end

function v = spval1(sp,x,left)
%SPVAL1 Evaluate univariate function in B-form.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);
%  If necessary, sort XS:
tosort = 0;
if any(diff(xs)< 0)
   tosort = 1; [xs,ix] = sort(xs);
end

%  Take apart spline:
[t,a,n,k,d] = spbrk(sp);
%  If there are no points to evaluate at, return empty matrix of appropriate
%  size:
if lx==0, v = zeros(d,0); return, end

%  Otherwise, augment the knot sequence so that first and last knot each
%  have multiplicity  >= K . (AUGKNT would not be suitable for this
%  since any change in T must be accompanied by a corresponding change
%  in A.)

index = find(diff(t)>0); addl = k-index(1); addr = index(length(index))-n;
if ( addl>0 | addr>0 )
   npk = n+k; t = t([ones(1,addl) 1:npk npk(ones(1,addr))]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr;
end

% For each data point, compute its knot interval:
if isempty(left)|left(1)~='l'
   index = max(sorted(t(1:n),xs),k);
else
   index = fliplr(max(n+1-sorted(-fliplr(t(k+1:n+k)),-xs),k));
end

% Now, all is ready for the evaluation.
if  k>1  % carry out in lockstep the first spline evaluation algorithm
         % (this requires the following initialization):
   dindex = reshape(repmat(index,d,1),d*lx,1);
   tx =reshape(t(repmat([2-k:k-1],d*lx,1)+repmat(dindex,1,2*(k-1))),d*lx,2*(k-1));
   tx = tx - repmat(reshape(repmat(xs,d,1),d*lx,1),1,2*(k-1));
   dindex = reshape(repmat(d*index,d,1)+repmat([1-d:0].',1,lx),d*lx,1);
   b = repmat([d*(1-k):d:0],d*lx,1)+repmat(dindex,1,k);
   a = a(:); b(:) = a(b);

   % (the following loop is taken from SPRPP)

   for r = 1:k-1
      for i = 1:k-r
         b(:,i) = (tx(:,i+k-1).*b(:,i)-tx(:,i+r-1).*b(:,i+1)) ./ ...
                  (tx(:,i+k-1)    -    tx(:,i+r-1));
      end
   end

   v = reshape(b(:,1),d,lx);
else     % the spline is piecewise constant, hence ...
   v = a(:,index);
end

if tosort>0, v(:,ix) = v; end

% Finally, zero out all values for points outside the basic interval:
index = find(x<t(1)|x>t(n+k));
if ~isempty(index)
   v(:,index) = zeros(d,length(index));
end
v = reshape(v,d*mx,nx);

function [knots,k] = aptknt(tau,k)
%APTKNT Acceptable knot sequence
%
%   APTKNT(TAU,K)  returns, for a given nondecreasing sequence TAU with
%   TAU(i) < TAU(i+K-1), all i, a knot sequence, KNOTS, for which the 
%   Schoenberg-Whitney conditions
%
%            KNOTS(i) <  TAU(i)  <  KNOTS(i+K) , i=1:length(TAU)
%
%   hold (with equality only for the first or last knot), ensuring that
%   the space of splines of order 
%              K  :=  min(K,length(TAU))  
%   with knot sequence KNOTS has a unique interpolant to arbitrary data
%   at the data sites TAU; the K used is, optionally, returned.
%
%   For example, for strictly increasing  x , and given corresponding  y ,
%
%      sp = spapi(aptknt(x,k),x,y);
%
%   gives a spline f  of order  min(k,length(x))  satisfying f(x(i)) = y(i),
%   all i (and the same result is obtained by spapi(k,x,y) ).
%   Be aware, though, of the fact that, for highly nonuniform  x , the 
%   determination of this spline can be ill-conditioned, leading possibly
%   to very strange behavior away from the interpolation points.
%
%   At present, the knot sequence chosen here is the initial guess used for
%   the iterative determination of the `optimal' knots in OPTKNT.
%
%   See also AUGKNT, AVEKNT, NEWKNT, OPTKNT.

%   Carl de Boor 14feb99, 07may00 (optionally return k as it may have changed)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2000/05/13 21:44:24 $

% If  tau(1) <= ... <= tau(n) with no more than  k-2  consecutive equalities,
% and  n>k , then the output  xi = aveknt(tau,k)  is strictly increasing and,
% for any a<tau(1), tau(n)<b, the output  knots = augknt([a xi b],k)  satisfies
% the above Schoenberg-Whitney conditions wrto  tau . 
%
% Indeed, then  
%                  knots(1:k) = a < tau(1) <= ... <= tau(k),
% while, for  i=1:n-k,
%           knots(k+i) =  xi(i) = (tau(i+1)+...+tau(i+k-1))/(k-1), 
% hence (using the fact that at most k-1 consecutive tau's can be equal)
%                 tau(i) < knots(k+i) < tau(i+k) ,   i=1:n-k ,
% and, finally,
%               tau(n-k+1) <= ... <= tau(n) < b = knots(n+[1:k]).  
% Letting now  a -->  tau(1)  and  b --> tau(end)  will not change any of these
% inequalities, except those involving the first and last data site may not
% be strict any more. But that is ok since these will be the endpoints of the
% corresponding basic interval, hence only right, respectively, left limits
% matter there.

n = length(tau); 
if n<2, error('There must be at least two (distinct) sites.'), end

k = max(1,min(k,n)); dtau = diff(tau);
if any(dtau<0)
   error('The site sequence must be nondecreasing.'), end
if k==1 % simply use midpoints between data sites
   if ~all(dtau)
      error('For k==1, the site sequence must be strictly increasing.')
   end
   knots = [tau(1) tau(1:n-1)+dtau/2 tau(n)];
else
   if any(tau(k:n)==tau(1:n-k+1))
      error(sprintf('No more than %g consecutive site(s) may coincide.',k-1))
   end
   knots = augknt([tau(1) aveknt(tau,k) tau(end)],k);
end

function pointer = sorted(meshsites, sites)
%SORTED Locate sites with respect to meshsites.
%
%   Given sequences MESHSITES and SITES, returns the index sequence
%   POINTER  (a  r o w  vector)  for which
%
%   POINTER(j) = #{ i : MESHSITES(i)  <=  sort(SITES)(j) },  all  j .
%
%   Thus, if both MESHSITES and SITES are nondecreasing, then
%
%   MESHSITES(POINTER(j))  <= SITES(j) < MESHSITES(POINTER(j)+1) ,
%
%   with the obvious interpretations when  POINTER(j) = 0  or
%   POINTER(j)+1 = length(MESHSITES) + 1 .
%
%   For example, 
%
%       sorted( [1 2 3 4], [0 1 2.1 2.99 3.5 4 5])
%
%   specifies 1:4 as MESHSITES and [0 1 2.1 2.99 3.5 4 5] as SITES and
%   gives the output  [0 1 2 2 3 4 4] .
%
%   See also PPUAL, SPVAL.

%   C. de Boor latest change: May 27, 1989
%   cb : 20 oc 94 : correct help
%   cb :  9 may '95 (use .' instead of ')
%   cb : 19 nov 95 (correct help; bring back to original, simple form since
%                   SORT doesn't reorder like terms any more)
%   cb : 28aug99 (point --> site, update help)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.9 $

[ignored,index] = sort([meshsites(:).' sites(:).']);
pointer = find(index>length(meshsites))-[1:length(sites)];

function [nbo,rows,ncols,last,blocks] = bkbrk(blokmat)
%BKBRK Part(s) of an almost block-diagonal matrix.
%
%   [NB,ROWS,NCOLS,LAST,BLOCKS] = BKBRK(BLOKMAT)
%
%   returns the details of the almost block diagonal matrix contained
%   in BLOKMAT, with ROWS and LAST NB-vectors, and BLOCKS a matrix
%   of size  SUM(ROWS)-by-NCOLS.
%
%   BKBRK(BLOKMAT)
%
%   returns nothing, but all parts are printed.
%
%   See also SPCOL, SLVBLK.

%   Carl de Boor 28 mar 90
%   cb :  1 feb 91: correct misspelling
%   cb : 16 jan 96: print only when no output arguments
%   cb : 17 mar 96: print out correct number of columns
%   cb : 5 oct 97 (replace use of ANS, to help compilation)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.10 $

if blokmat(1)==41 % data type number for the spline block format is 41
   % Here are the details of this particular sparse format:
   % The matrix is sum(ROWS)-by-sum(LAST).
   % There are NB blocks. The i-th block has ROWS(i) rows and NCOLS columns.
   % The first column of the (i+1)st block is exactly LAST(i) columns to the
   % right of the first column of the i-th block.
   nb = blokmat(2);
   rows = blokmat(2+[1:nb]);
   ncols = blokmat(2+nb+1);
   last = blokmat(3+nb+[1:nb]);
   blocks = reshape(blokmat(3+2*nb+[1:sum(rows)*ncols]),sum(rows),ncols);

elseif blokmat(1)==40 % data type number for general almost block diagonal
                      % format is 40;
   nb = blokmat(2);
   rows = blokmat(2+[1:nb]);
   cols = blokmat(2+nb+[1:nb]);
   last = blokmat(2+2*nb+[1:nb]);
   row = cumsum([0,rows]);
   ne = sum(rows);ncols = max(cols);
   len = rows.*cols;
   index = cumsum([2+3*nb len]);
   blocks = zeros(ne,ncols);
   for j=1:nb
      block = reshape(blokmat(index(j)+[1:len(j)]),rows(j),cols(j));
      blocks(row(j)+[1:row(j+1)],[1:cols(j)]) = block;
   end
else
   error('The argument does not appear to be an almost block diagonal matrix.')
end

if nargout==0 % print out the blocks
   if blokmat(1)==41 % generate COLS
      temp = cumsum([0 last]); temp = temp(nb+1)-temp;
      cols = min([temp(1:nb);ncols(1,ones(1,nb))]);
   end
   rowsum = cumsum([0 rows]);
   for j=1:nb
      fprintf(['block ',int2str(j),' has ',int2str(rows(j)),' row(s)\n'])
      disp(blocks(rowsum(j)+[1:rows(j)],1:cols(j)))
    fprintf([' next block is shifted over ',int2str(last(j)),' column(s)\n\n'])
   end
else
   nbo = nb;
end

function [augknot,addl] = augknt(knots,k,mults)
%AUGKNT Augment a knot sequence.
%
%   AUGKNT(KNOTS,K) returns a nondecreasing and augmented knot 
%   sequence which has the first and last knot with exact multiplicity K .  
%   (This may actually shorten the knot sequence.)  
%
%   [AUGKNOT,ADDL] = AUGKNT(KNOTS,K) also returns the number of knots
%   added on the left.  (This may be negative.)
%
%   AUGKNOT = AUGKNT(KNOTS,K,MULTS) returns an augmented knot sequence
%   that, in addition, contains each interior knot MULTS times.  If
%   MULTS has exactly as many entries as there are interior knots,
%   then the j-th one (in the ordered sequence) will be repeated MULTS(j)
%   times.  Otherwise, the uniform multiplicity MULTS(1) is used, whose
%   default value is 1 .  If the sequence of interior knots in KNOTS is
%   *strictly* increasing, then this ensures that the splines of order K
%   with knot sequence AUGKNOT satisfy K-MULTS(j) smoothness conditions
%   across the j-th interior break, all j .
%
%   For example, the statement
%
%      ch = spapi(augknt(x,4,2), [x x], [y dy]);
%
%   constructs the piecewise cubic Hermite interpolant, i.e., the function  s
%   in  ch  satisfies  s(x(i)) = y(i),  (Ds)(x(i)) = dy(i), all i (assuming 
%   that  x  is strictly increasing).
%
%   See also SPALLDEM.

%   Carl de Boor 25 feb 89
%   cb:  5 apr 95 (individual multiplicities)
%   cb: 17 nov 95 (permit unordered input, and add additional input checking)
%   cb: 22sep96 (use BRK2KNT)
%   cb: 06oct97, 04may98 (improve and standardize the help)
%   cb: 11dec98 (use repmat)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.11 $  $Date: 2000/04/25 14:13:52 $

if nargin<3
   if (length(k)>1|k<1)
      error('The second argument should be a single natural number.'), end
   mults = 1;
end

dk = diff(knots);
if ~isempty(find(dk<0)), knots = sort(knots); dk = diff(knots); end

augknot=[];
j=find(dk>0); if isempty(j)
   error('The knot sequence should contain more than one point.'), end
addl = k-j(1);

interior = (j(1)+1):j(end);
%   % make sure there is a multiplicity assigned to each interior knot:
if length(mults)~=length(interior), mults = repmat(mults(1),size(interior)); end

augknot = brk2knt(knots([1 interior end]), [k mults k]);

function fprime = fnder(f,dorder)
%FNDER Differentiate a function.
%
%   FNDER(F)
%
%   returns the (representation of the) first derivative of the
%   univariate function contained in F (and in the same form).  
%
%   FNDER(F,DORDER) returns the DORDER-th derivative, with DORDER expected
%   to be of the form [d1,...,dm] in case the function in F is m-variate,
%   and, for each i=1,..,m,  di  an integer to indicate that the function
%   in F is to be differentiated di-fold with respect to its i-th argument.
%   Here, di may be negative, resulting in di-fold integration with respect
%   to the i-th argument.
%
%   This command does not work for rational splines; for them, use FNTLR 
%   instead.
%
%   For example,
%
%   fnval( fnder( sp, 2), 3.14 );
% 
%   gives the value at 3.14 of the function in sp, while
%
%   sp0 = fnint( fnder( sp ) );
%
%   gives a function that differs from sp by some constant only (namely, by
%   its value at 0).
%
%   See also FNDIR, FNTLR, FNINT.

%   C. de Boor latest change: Feb.25, 1989
%   cb 20 oct 94 : remove extraneous argument from call to  spmak
%   cb  4 nov 94 : include differentiation of bivariate tensor product splines
%                  in both B-form and ppform.
%   cb : 9 may '95 (use .' instead of '), 29aug96 (include BB-form)
%   cb : 5 oct 97 (replace use of ANS, to help compilation)
%   cb : 26oct97 (treat also multivariate functions correctly;
%                 permit DORDER to be negative to do integration)
%   cb : 05may98 (standardize help; accept earlier versions of forms)
%   cb : 11dec98 (replicate via repmat)
%   cb : 23jan00 (refer to FNTLR for rational splines)
%   cb : 09apr00 (deal with fact that matlab drops trailing singleton dims)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.11 $

if nargin<2, dorder=1; end

if ~isstruct(f), f = fn2fm(f); end

switch f.form(1)
case 'p' % the function is in ppform:
   [breaks,coefs,l,k,d]=ppbrk(f);
   if iscell(breaks) % the function is multivariate
      m = length(k);
      if length(dorder)~=m
         error(['DORDER should be a ' num2str(m) '-vector.']), end
      sizec = [d,l.*k]; %size(coefs);
      for i=m:-1:1
         dd = prod(sizec(1:m));
         dpp = fnderp(ppmak(breaks{i},reshape(coefs,dd*l(i),k(i)),dd), ...
                      dorder(i));
         breaks{i} = dpp.breaks; sizec(m+1) = dpp.pieces*dpp.order;
         coefs = reshape(dpp.coefs,sizec);
         if m>1
             coefs = permute(coefs,[1,m+1,2:m]);
             sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = ppmak(breaks,coefs,sizec);
   else
      fprime = fnderp(f,dorder);
   end

case 'B' % the function is in B-form or BB-form;
                           % omit trivial B-spline terms.
   [knots,coefs,n,k,d]=spbrk(f);
   if iscell(knots)       % the function is multivariate
      m = length(knots);
      if length(dorder)~=m
         error(['DORDER should be a ' num2str(m) '-vector.']), end
      sizec = [d,n];% size(coefs);
      for i=m:-1:1
         dsp = fnderb(spmak(knots{i},...
            reshape(coefs,prod(sizec(1:m)),sizec(m+1))),dorder(i));
         knots{i} = dsp.knots; sizec(m+1) = dsp.number; 
         coefs = reshape(dsp.coefs,sizec); 
         if m>1
            coefs = permute(coefs,[1,m+1,2:m]);
            sizec(2:m+1) = sizec([m+1,2:m]);
         end
      end
      fprime = spmak(knots,coefs,sizec);
   else
      fprime = fnderb(f,dorder);
   end
case 'r'
   error('This command does not work for rational splines. Use FNTLR instead.')
otherwise
   error('F does not appear to describe a function.')
end

function fprime = fnderp(f,dorder)
%FNDERP Differentiate a univariate function in ppform.
[breaks,coefs,l,k,d]=ppbrk(f);
if k<=dorder
   fprime=ppmak([breaks(1) breaks(l+1)],zeros(d,1));
elseif dorder<0    % we are to integrate
   fprime = f;
   for j=1:(-dorder)
      fprime = fnint(fprime);
   end
else
   knew=k-dorder;
   for j=k-1:-1:knew
      coefs=coefs.*repmat([j:-1:j-k+1],d*l,1);
   end
   fprime=ppmak(breaks,coefs(:,1:knew),d);
end

function fprime = fnderb(f,dorder)
%FNDERB Differentiate a univariate function in B-form.

[t,a,n,k,d]=spbrk(f);
if k<=dorder
   fprime=spmak(t,zeros(d,n));
elseif dorder<0    % we are to integrate
   fprime = f;
   for j=1:(-dorder)
      fprime = fnint(fprime);
   end
else
   knew=k-dorder;
   for j=k-1:-1:knew
      tt=t(j+1+[0:n])-t(1:n+1); z=find(tt>0); nn=length(z);
      temp=(diff([zeros(1,d);a.'; zeros(1,d)])).';
      a=temp(:,z)./repmat(tt(z)/j,d,1);
      t=[t(z) t(n+2:n+j+1)]; n=nn;
   end
   fprime=spmak(t,a);
end

function t = brk2knt(breaks,mults)
%BRK2KNT Breaks with multiplicities into knots.
%
%   T = BRK2KNT(BREAKS,MULTS)
%   returns the sequence T in which, for each i, BREAKS(i) is repeated
%   MULTS(i) times. If BREAKS is strictly increasing, then T is the
%   knot sequence in which each BREAKS(i) occurs exactly MULTS(i) times.
%
%   If MULTS is to be constant, then only that constant value need
%   be given.
%
%   [T,INDEX] = BRK2KNT(BREAKS,MULTS)
%   also returns INDEX = [1 find(diff(T)>0)-1] . If all multiplicities are
%   positive, then, for all j, INDEX(j) is the first place in T at which
%   BREAKS(j) appears.
%
%   For example,
%
%      t = [1 1 2 2 2 3 4 5 5];  [xi,m] = knt2brk(t);  tt = brk2knt(xi,m);
%
%   gives  [1 2 3 4 5]  for xi , [2 3 1 1 2] for m , and t for tt.
%
%   See also KNT2BRK, KNT2MLT, AUGKNT.

%   cb 26aug96
%   cb 06oct97 (improve the help)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.7 $

s = sum(mults);
if s==0
   t = [];
else
   li = length(breaks);
      % make sure there is a multiplicity assigned to each break,
      % and drop any break whose assigned multiplicity is not positive.
   if length(mults)~=li, mults = mults(ones(1,li));
   else
      fm = find(mults<=0);
      if ~isempty(fm), breaks(fm)=[]; mults(fm)=[]; li = length(breaks); end
   end
   mm = zeros(1,s);
   mm(cumsum([1 reshape(mults(1:li-1),1,li-1)])) = ones(1,li);
   t = breaks(cumsum(mm));
end

function [out1,coefs,n,k,d] = spbrk(sp,part)
%SPBRK Part(s) of a B-form or a BBform.
%
%   [KNOTS,COEFS,N,K,D] = SPBRK(SP) breaks the B-form in SP into its parts and 
%   returns as many of them as are specified by the output arguments.
%
%   OUT1 = SPBRK(SP,PART) returns the part specified by the string PART which 
%   may be (the beginning character(s) of) one of the following strings: 
%   'knots' or 't', 'coefs', 'number', 'order', 'dimension', 'interval' 
%
%   If PART is the 1-by-2 matrix [A,B], the restriction/extension of the spline
%   in SP to the interval with endpoints A and B is returned, in the same form.
%
%   SPBRK(SP) returns nothing, but prints out all the parts.
%
%   See also PPBRK, FNBRK, RSBRK, RPBRK.

%   Carl de Boor 25 feb 89
%   cb : 28 aug 94 (enlarge the options)
%   cb : 20 nov 95 (enlarge the options), 29aug96 (add BB-form)
%   cb : 18oct97 (change to structure, also cover multivariate spline)
%   cb : 22apr98 (accept array version of form; standardize help) 
%   cb : 09jan00 (minor  improvements in programming style and help)
%   cb : 12feb00 (check that second argin is string)
%   cb : 16apr00 (provide change of basic interval)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.13 $

if ~isstruct(sp)
  if sp(1)~=11&sp(1)~=12
     error('The input array does not seem to describe a function in B-form.')
  else
     di=sp(2);ni=sp(3);
     ci=reshape(sp(3+[1:di*ni]),di,ni);
     kk=sp(4+di*ni);ki=sp(4+di*ni+[1:kk+ni]);
     sp = spmak(ki,ci);
  end
end

if sp.form(1)~='B'
   error('The input does not seem to describe a spline in B-form.')
end
if nargin>1
   if nargout>1, error('Too many output arguments for the given input.')
   elseif ~isstr(part)
      if iscell(part)  % we must be dealing with a tensor-product spline
         c = sp.coefs; knots = sp.knots; m = length(knots);
         sizec = size(c);
         if length(sizec)~=m+1 % in trouble because of trailing singleton dims
            sizec = [sp.dim,sp.number]; c = reshape(c,sizec);
         end
         for i=m:-1:1
            dd = prod(sizec(1:m));
            spi = spcut(spmak(knots{i},reshape(c,dd,sp.number(i))), part{i});
            knots{i} = spi.knots; sizec(m+1) = spi.number;
            c = reshape(spi.coefs,sizec);
            if m>1
               c = permute(c,[1,m+1,2:m]);
               sizec(2:m+1) = sizec([m+1,2:m]);
            end
         end
         out1 = spmak(knots,c,sizec);

      else             % we must be dealing with a univariate spline
         out1 = spcut(sp,part);
      end
   else
      switch part(1)
      case 'd',       out1 = sp.dim;
      case 'n',       out1 = sp.number;
      case {'k','t'}, out1 = sp.knots;
      case 'o',       out1 = sp.order;
      case 'c',       out1 = sp.coefs;
      case 'i', % this must be treated differently in multivariate case
         if iscell(sp.knots)
            for i=length(sp.knots):-1:1  % loop backward to avoid redef.
               out1{i} = sp.knots{i}([1 end]);
            end
         else
           out1 = sp.knots([1 end]);
         end
      otherwise
         error(['''',part,''' is not part of a B-form.'])
      end
      return
   end
else
   if nargout==0
     if iscell(sp.knots) % we have a multivariate spline and, at present,
                         % I can't think of anything clever to do; so...
       disp(sp)
     else
       disp('knots(1:n+k)'),disp(sp.knots),
       disp('coefficients(d,n)'),disp(sp.coefs),
       disp('number n of coefficients'),disp(sp.number),
       disp('order k'),disp(sp.order),
       disp('dimension d of target'),disp(sp.dim),
     end
   else
    out1 = sp.knots; coefs = sp.coefs; n = sp.number; k = sp.order; d = sp.dim;
   end
end
function out1 = spcut(sp,interv)
%SPCUT change the basic interval
   if isempty(interv), out1 = sp; return, end
   sizei = size(interv);
   if sizei(2)>1 % we are to change the basic interval
      tl = interv(1,1); tr = interv(1,2);
      if tl==tr
         warning('No changes made since the given end points are equal.')
         out1 = sp; return
      end
      if tl>tr, tl = tr; tr = interv(1); end

      index = sorted(sp.knots,[tl,tr]); mults = knt2mlt(sp.knots);
      if tl<sp.knots(1),      m1 = 1;
      elseif tl==sp.knots(1), m1 = 0;
      else,                   m1 = sp.order;                    
         if tl==sp.knots(index(1)), m1 = m1-mults(index(1))-1; end
      end
      if tr>sp.knots(end),      m2 = 1;
      elseif tr==sp.knots(end), m2 = 0;
      else,                     m2 = sp.order;                    
         if tr==sp.knots(index(2)), m2 = m2-mults(index(2))-1; end
      end
      sp = fnrfn(sp, [repmat(tl,1,m1),repmat(tr,1,m2)]);
      index = sorted(sp.knots,[tl tr]);
      if sp.knots(end)>tr
         sp = spmak(sp.knots(1:index(2)),sp.coefs(:,1:(index(2)-sp.order)));
      end
      if sp.knots(1)<tl
         sp = spmak(sp.knots(index(1)-sp.order+1:end), ...
                    sp.coefs(:,index(1)-sp.order+1:end));
      end
      out1 = sp;
   else 
      error('The second argument should be a string or an interval.')
   end

   function spline = spmak(knots,coefs,sizec)
%SPMAK Put together a spline in B-form.
%
%   SPMAK(KNOTS,COEFS) puts together a spline from the knots and 
%   coefficients. COEFS should be d-by-n , with  d  the number of
%   components in each of the  n  spline coefficient supplied.  
%   The order of the spline is inferred as  k := length(KNOTS) - n .
%   Knot multiplicity is held to <= k , with the coefficients
%   corresponding to a B-spline with trivial support ignored.
%
%   SPMAK  will prompt for KNOTS and COEFS.
%
%   If KNOTS is a cell array of length  m , then COEFS is expected to
%   be an m- or (m+1)-dimensional array, and the corresponding m-variate
%   d-valued tensor product spline is returned, with  d  the first 
%   dimension of COEFS in case COEFS is (m+1)-dimensional, and 1 otherwise.
%   FEWER CHECKS ARE CARRIED OUT in this case.  
%
%   SPMAK(KNOTS,COEFS,SIXEC) uses SIZEC to specify the intended array
%   dimensions of COEFS, and may be needed for proper interpretation
%   of COEFS in case one or more of its trailing dimensions is a singleton
%   and thus COEFS appears to be of lower dimension.
%
%   For example, if the intent is to construct a 2-vector-valued bivariate
%   polynomial on the rectangle [-1 .. 1] x [0 .. 1], linear in the first
%   variable and constant in the second, say
%      coefs = zeros([2 2 1]); coefs(:,:,1) = [1 0;0 1];
%   then the straightforward
%      sp = spmak({[-1 -1 1 1],[0 1]},coefs);
%   will fail, producing a scalar-valued function of (illegal) order [2 0],
%   as will
%      sp = spmak({[-1 -1 1 1],[0 1]},coefs,size(coefs));
%   while proper use of that third argument, as in
%      sp = spmak({[-1 -1 1 1],[0 1]},coefs,[2 2 1]);
%   will succeed.
%
%   See also SPBRK, RSMAK, PPMAK, RPMAK, FNBRK.

%   Carl de Boor 9 dec 89
%   cb : February 4, 1991 (disallow the check for d>n)
%   cb : 20 oct 1994 (enforce knot multiplicity <= k )
%   cb :  9 may '95 (use .' instead of ')
%   cb : 15aug96 (check that basic interval is nontrivial)
%   cb : 18oct97 (change to structure, also cover multivariate spline)
%   cb : 07dec97 (reorder form for backward consistency, ensure knots is 1-row)
%   cb : 10may98 (standardize the help)
%   cb : 10apr00 (introduce optional third argument for multidimensional case)
%   cb : 23jun00 (improve the help)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.10 $

if nargin==0;
   knots = input('Give the vector of knots  >');
   coefs = input('Give the array of B-spline coefficients  >');
end
if iscell(knots) % we are putting together a tensor-product spline

   m = length(knots); if nargin<3, sizec = size(coefs); end
   switch length(sizec)
      case m, sizec = [1, sizec]; coefs = reshape(coefs,sizec);
      case m+1
      otherwise
         error('COEFS must be a ([1+]length(KNOTS))-dimensional array.')
   end
   d = sizec(1); n = sizec(2:end);
   for j=m:-1:1
      k(j) = length(knots{j})-n(j);
      if k(j)<=0, error('There should be more knots than coefficients.'), end
      knots{j} = reshape(knots{j},1,k(j)+n(j));
   end

else            % we are putting together a univariate spline

   [d,n] = size(coefs);
   if isempty(coefs), error('The coefficient sequence is empty.'), end
   
   k = length(knots)-n;
   if k<=0, error('There should be more knots than coefficients.'), end
   diffk = diff(knots);
   if any(diffk<0)
      error('The knot sequence should be nondecreasing.')
   end
   if knots(1)==knots(end)
      error('The extreme knots should be different.')
   end
   
   index = find(knots(k+[1:n])-knots(1:n)>0); % throw out trivial B-splines
   if length(index)<n
      coefs = coefs(:,index); knots = knots([index n+[1:k]]); n = length(index);
   end
   knots = reshape(knots,1,n+k);
end

spline.form = 'B-';
spline.knots = knots;
spline.coefs = coefs;
spline.number = n;
spline.order = k;
spline.dim = d;
% spline = [11 d n coefs(:).' k knots(:).'];

function tstar = aveknt(t,k)
%AVEKNT Knot averages.
%
%   AVEKNT(T,K)  returns the averages of successive  K-1  knots, i.e., 
%   the points 
%
%       TSTAR(i) = ( T_{i+1} + ... + T_{i+K-1} ) / (K-1)
%
%   recommended as good interpolation point choices when interpolating
%   from  S_{K,T} .
%
%   For example, with  k  and the increasing sequence  breaks  given,
%   the statements
%
%      t = augknt(breaks,k); x = aveknt(t);
%      sp = spapi( t , x, sin(x) );
%
%   provide a spline interpolant to the sine function on the interval
%   [breaks(1) .. breaks(end)] .
%
%   See also SPAPIDEM, OPTKNT, APTKNT, CHBPNT.

%   Carl de Boor 24 may 89
%   cb : 22 mar 91 (program around a MATLAB discontinuity)
%   cb :  9 may 95 (use .' instead of ')
%   cb :  5 mar 96 (use reshape etc.)
%   cb :  5 oct 97 (replace use of ANS, to help compilation)
%   cb :  4 may 98 (standardize help)
%   cb : 11 dec 98 (use repmat)
%   cb : 14 feb 99 (handle tau of length k)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc.
%   $Revision: 1.11 $

t = t(:); n = length(t)-k;
if k<2, error('The second argument must be at least 2 .')
elseif n<0, error('There must be at least K knots.')
elseif k==2, tstar = reshape(t(1+[1:n]),1,n);
else
   temp = repmat(t,1,k-1);
   temp = sum(reshape([temp(:);zeros(k-1,1)],n+k+1,k-1).')/(k-1);
   tstar = temp(1+[1:n]);
end

function x = slvblk(blokmat,b,w)
%SLVBLK Solve almost block-diagonal linear system.
%
%   SLVBLK(BLOKMAT,B)  returns the solution (if any) of the linear system  
%   A*X=B, with the matrix A stored in BLOKMAT in the spline almost block
%   diagonal form (as generated, e.g., in SPCOL).
%
%   If the system is overdetermined (i.e., has more equations than
%   unknowns), the least-squares solution is returned.  In that case, it
%   may be useful to use the optional third argument, i.e., to use the
%   following:
%
%   SLVBLK(BLOKMAT,B,W)  returns the vector X that minimizes the 
%   w e i g h t e d  l_2 sum
%
%      sum_j W(j)*( (A*X-B)(j) )^2 .
%
%   The default for W is the sequence [1,1,1,...].
%
%   See also SPCOL, SPAPS, SPAPI, SPAP2.

%   Carl de Boor 14 jan 90
%   cb : Nov.29, 1990 (correct signum problem)
%   cb : Feb. 1, 1991 (correct signum problem further)
%   cb : 21 Jun, 1992 (change .. to full ellipsis ...)
%   cb : 15 sep, 1994 (add optional weight)
%   cb : 11 nov, 1994 (remove wasteful multiplication by identity)
%   cb : 13sep96 (add sparse matrix option)
%   cb : 10may98 (standardize the help)
%   cb : 13jul98 (initialize ELIM to avoid confusion with any ELIM.M)
%   Copyright 1987-2000 by C. de Boor and The MathWorks, Inc. 
%   $Revision: 1.9 $

% If BLOKMAT is sparse, handle the problem sparsely:
if issparse(blokmat)
   if nargin>2
      n = length(w); spw = sparse(1:n,1:n,sqrt(w));
      x = (spw*blokmat)\(spw*b);
   else
      x = blokmat\b;
   end
   return
end

% get the basic information
[nb,rows,ncols,last,blocks] = bkbrk(blokmat);

ne = sum(rows);nu = sum(last);
if any(cumsum(rows)<cumsum(last))|any(last>ncols)
   error('The coefficient matrix has a nontrivial nullspace.')
end

[brow,bcol] = size(b);
if(ne~=brow)
   error('Matrix and right side are incompatible.')
end

blocks = [blocks b];
ccols = ncols+bcol;
if nargin>2, w = sqrt(w); blocks = (w(:)*ones(1,ccols)).*blocks; end

f = 1; l = 0; elim = 0;
for j=1:nb
   if (f<=l) % shift the rows still remaining from previous block
      blocks(f:l,:) = ...
         [blocks(f:l,elim+1:ncols) zeros(l+1-f,elim),blocks(f:l,ncols+1:ccols)];
   end
   l = l+rows(j);

   elim = last(j);
   % ideally, one would now use
   %   [q,r] = qr(blocks(f:l,1:elim));
   % followed up by
   %   blocks(f:l,:) = q'*blocks(f:l,:);
   %   f = f+elim;
   % but, unfortunately, this generates the possibly very large square matrix q
   % The unhappy alternative is to do the elimination explicitly here, using
   % Householder reflections (and an additional inner loop):
   for k=1:elim
      a = norm(blocks(f:l,k));
      vv = abs(blocks(f,k))+a;
      c = vv*a;
      if blocks(f,k)<0, vv = -vv; end
      q = [vv;blocks(f+1:l,k)];
      blocks(f:l,:) = ...
       blocks(f:l,:)-((q/c)*ones(1,ccols)).*(ones(l+1-f,1)*(q'*blocks(f:l,:)));
      f = f+1;
   end
end

% now we are ready for back-substitution
x = zeros(f-elim-1+ncols,bcol);

for j=nb:-1:1
   elim = last(j); l = f-1; f = f-elim;
   % here is another occasion where empty matrices of various sizes would help;
   % instead, use an if statement:
   if elim<ncols, blocks(f:l,ncols+1:ccols) = blocks(f:l,ncols+1:ccols) ...
                    - blocks(f:l,elim+1:ncols)*x(f-1+[elim+1:ncols],:); end
   x(f:l,:) = blocks(f:l,1:elim) \ blocks(f:l,ncols+1:ccols);
end
x = x(1:nu,:);
