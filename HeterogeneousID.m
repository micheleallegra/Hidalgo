%mex FFLAGS="-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -ffree-line-length-0"  DPA.F90 rnkpar.f90

function out=AdvancedDensityPeaks(varargin)

K=2;
q=3;
run_K=1;
zeta = 0.8;
Niter= 10000;
Nreplica = 1;
type = 1;
ass_th=0.8;

if(nargin==0)
  error('At least input data required.');
end

if(nargout==0)
  error('At least one output required.');
end

ctr=0;

if(isnumeric(varargin{ctr+1}))
    X=varargin{ctr+1};
    N = size(X,1);
    ctr=ctr+1;
elseif(ischar(varargin{ctr+1}))
    if(varargin{ctr+1}=='Distance')
       type=0;
       if(size(varargin{ctr+2},1)==size(varargin{ctr+2},2))
         D=varargin{ctr+2};
         N = size(D,1);
         mask = tril(true(size(D)),-1);
         D=D(mask);
         ctr=ctr+2;
       else
         ss=['input #', num2str(ctr+2), ' must be a square distance matrix'];
         error(ss);
       end
    else
      ss=['Set input #', num2str(ctr+1), ' to the input coordinate matrix, or set input#', num2str(ctr+1),'to ''Distance'' and input', num2str(ctr+2),'to the input distance matrix'];
      error(ss);
    end
else
   ss=['Set input #', num2str(ctr+1), ' to the input coordinate matrix, or set input#', num2str(ctr+1),'to ''Distance'' and input', num2str(ctr+2),'to the input distance matrix'];
   error(ss);
end


while(ctr < nargin)

  option=varargin{ctr+1};

  switch(option)

    case ('K')
      if(nargin < ctr+2)
        ss=['Specify k to use as argument #',num2str(ctr+2)];
        error(ss);
      elseif(isnumeric(varargin{ctr+2})==0)
        ss=['Input #',num2str(ctr+2),'must be an integer number'];
        error(ss);
      else
        sk = length(varargin{ctr+2});
        if(sk > 1)
          ss=['Input #',num2str(ctr+2),'must be an integer number'];
          error(ss);
        end
        K=varargin{ctr+2};
	run_K=0;
      end
      ctr=ctr+2;

  case ('q')
      if(nargin < ctr+2)
        ss=['Specify q to use as argument #',num2str(ctr+2)];
        error(ss);
      elseif(isnumeric(varargin{ctr+2})==0)
        ss=['Input #',num2str(ctr+2),'must be an integer number'];
        error(ss);
      else
        sq = length(varargin{ctr+2});
        if(sq > 1)
          ss=['Input #',num2str(ctr+2),'must be an integer number'];
          error(ss);
        end
        q=varargin{ctr+2};
      end
      ctr=ctr+2;


  case ('zeta')
      if(nargin < ctr+2)
        ss=['Specify q to use as argument #',num2str(ctr+2)];
        error(ss);
      elseif(isnumeric(varargin{ctr+2})==0)
        ss=['Input #',num2str(ctr+2),'must be an integer number'];
        error(ss);
      else
        sz = length(varargin{ctr+2});
        if(sz > 1)
          ss=['Input #',num2str(ctr+2),'must be an integer number'];
          error(ss);
        end
        zeta=varargin{ctr+2};
      end
      ctr=ctr+2;

  case ('Niter')
      if(nargin < ctr+2)
        ss=['Specify number of iterations as argument #',num2str(ctr+2)];
        error(ss);
      elseif(isnumeric(varargin{ctr+2})==0)
        ss=['Input #',num2str(ctr+2),'must be an integer number'];
        error(ss);
      else
        snit = length(varargin{ctr+2});
        if(snit > 1)
          ss=['Input #',num2str(ctr+2),'must be an integer number'];
          error(ss);
        end
        Niter=varargin{ctr+2};
      end
      ctr=ctr+2;

  case ('Nreplicas')
      if(nargin < ctr+2)
        ss=['Specify numer of replicas as argument #',num2str(ctr+2)];
        error(ss);
      elseif(isnumeric(varargin{ctr+2})==0)
        ss=['Input #',num2str(ctr+2),'must be an integer number'];
        error(ss);
      else
        snr = length(varargin{ctr+2});
        if(snr > 1)
          ss=['Input #',num2str(ctr+2),'must be an integer number'];
          error(ss);
        end
        Nreplica=varargin{ctr+2};
      end
      ctr=ctr+2;
	    	
 otherwise
      ss=['unrecognized input ',option];
      error(ss);
 end
end

Lk=-1E10;
if(run_K)
	for K=1:8
		if(type==1)
			[d,p,Pi,L]=hidalgo(1,X',K,q,zeta,Niter,Nreplica);
		else
			[d,p,Pi,L]=hidalgo(0,D,K,q,zeta,Niter,Nreplica);
		end

		Pi = reshape(Pi,[N*K Nreplica]);
		Pi = reshape(Pi,[K N Nreplica]);
		[Lmax,rmax] = max(L(2,:));
		if(Lmax>Lk)
			Lk=Lmax;	
			out.d = d(:,rmax);
			out.p = p(:,rmax);
			out.Pi = Pi(:,:,rmax);
			[Pimax, Z] = max(Pi);
			unassigned = find(Pimax<ass_th);
			Z(unassigned)=0;
			out.Z = Z;
		else
			break;
		end
	end
	
else
	if(type==1)
		[d,p,Pi,L]=hidalgo(1,X',K,q,zeta,Niter,Nreplica);
	else
		[d,p,Pi,L]=hidalgo(0,D,K,q,zeta,Niter,Nreplica);
	end

	Pi = reshape(Pi,[N*K Nreplica]);
	Pi = reshape(Pi,[K N Nreplica]);

	[Lmax,rmax] = max(L(2,:));

	out.d = d(:,rmax);
	out.p = p(:,rmax);
	out.Pi = Pi(:,:,rmax);
        [Pimax, Z] = max(Pi);
        unassigned = find(Pimax<ass_th);
        Z(unassigned)=0;
	out.Z = Z;
        out.L = L(2,rmax);
end


end
