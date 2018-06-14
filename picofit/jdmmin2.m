
% [retfitinfo,yfit,amp]=jdmmin(DataX, DataY, fitinfo)
% minimization function, for TCSPC data, loosely based on simplex
% minimization principles.

% yfit is a matrix and amp is a vector such that yfit*amp is the
% overall fit function, while yfit(:,N) is the N-th function, and amp(N) is the N-th
% amplitude.

% [retfitinfo,yfit,amp]=kwwfit(DataX,DataY,fitinfo) ;
% [sqerror,yfit,AMP]=kwwfit_fun(fitinfo);

% approach is to evaluate the function at the edge and center values of an N-dimensional hypercube
% where N is the number of nonlinear parameters. 

global IRFx IRFy 

clear parmidx parmstepsz 

% in function version, will specify these in input instead of here
betastepsz=0.0000002;
t0stepsz=0.02; 
taustepsz=0.0002;

% part 1: determine how many exp and kww functions in fitinfo
numexp=0;
numkww=0;
%parmnum=1; % for keeping track of (0 corresponds to t0 parameter)

parmidx{1}=[0 0]; % corresponds to t0 shift parameter
parmstepsz(1)=t0stepsz;

for cnt=1:length(fitinfo.fun),
    if fitinfo.fun{cnt}.function=='exp'
        parmidx{end+1}=[cnt 1];
        parmstepsz(end+1)=taustepsz;
        numexp=numexp+1;
    end
    if fitinfo.fun{cnt}.function=='kww'
        parmidx{end+1}=[cnt 1]; % tau parameter
        parmstepsz(end+1)=taustepsz;
        parmidx{end+1}=[cnt 2]; % beta parameter
        parmstepsz(end+1)=betastepsz;
        numkww=numkww+1;
    end
end

% note: parmidx{1..N} containds the locations (indices) of all of the nonlinear parameters

% part 1.1: determine the number of nonlinear variables and the size of the parameter spzce (number of parameter sets)
numvar=1+numexp+2*numkww;
numeval=3^numvar;

% part 2: generate a list of potential moves corresponding to each parameter value to evaluate.
% this may involve using integer modulo arithmetic, giving -1,0,1 for each move value. Or alternatively,
% number each nonlinear parameter value 1..3^N
% Uses dec2base to generate base-3 numbers
moves=zeros([numeval numvar]);
for cnt=1:numeval
    movstr=dec2base(cnt-1,3,numvar);
    % movstr=movstr(end:-1:1); % for debugging?
    for cnt2=1:numvar
        moves(cnt,cnt2)=str2num(movstr(cnt2))-1;
    end
    if sum(moves(cnt,:).*moves(cnt,:))==0
        centerpoint=cnt; % center point of hypercube
    end
end

% save fitinfo, so we can compare to refined fit
fitinfo0=fitinfo;
% save square error
[sqerror,yfit,amp]=kwwfit_fun(fitinfo);
sqerr0=sqerror;

starttime=cputime;

%  part 3 main loop - searching for minimum
while 1

    % part 3.1, loop over set of parameters
    sqerrs=zeros([numeval 1]);
    for cnt=1:numeval
        % part 3.2: inner loop: loop over parameters. evaluate for each parameter set
        % note: make sure not to exceed bounds
        curparmset=moves(cnt,:);
        fittmp=fitinfo;
        for curvar=1:numvar
            % make changes to each "curval"
            if curvar==1 % for t0 parameter
                fittmp.t0curval=fittmp.t0curval+curparmset(1)*parmstepsz(1);
            else
                parmtmp=fittmp.fun{parmidx{curvar}(1)}.curval(parmidx{curvar}(2));
                fittmp.fun{parmidx{curvar}(1)}.curval(parmidx{curvar}(2))=parmtmp+curparmset(curvar)*parmstepsz(curvar);
            end
        end
        % evaluate and store in sqerrs(cnt)
        [sqerror,yfit,amp]=kwwfit_fun(fittmp);
        sqerrs(cnt)=sqerror;
        fitinfoarray(cnt)=fittmp;
    end % end for
    
    % find minimum sqerr
    [val,idx]=min(sqerrs);
    % part 4: if min is at the origin , then reduce step size, or "break"
    if idx==centerpoint
        disp('minimum reached');
        break;
    end
    % part 4.1: move based on minimum value. Update curval in each function in fitinfo
    % note: make sure not to exceed bounds
    %disp('possible move found');
    newfitinfo=fitinfoarray(idx); % candidate fitinfo
    for cnt=1:length(newfitinfo.fun),
        if newfitinfo.fun{cnt}.function=='exp'
            if newfitinfo.fun{cnt}.curval>newfitinfo.fun{cnt}.ubound
                %disp('exp over bound');
                newfitinfo.fun{cnt}.curval=newfitinfo.fun{cnt}.ubound;
            end
            if newfitinfo.fun{cnt}.curval<newfitinfo.fun{cnt}.lbound
                %disp('exp under bound');
                newfitinfo.fun{cnt}.curval=newfitinfo.fun{cnt}.lbound;
            end
        end
        if fitinfo.fun{cnt}.function=='kww'
            if newfitinfo.fun{cnt}.curval(1)>newfitinfo.fun{cnt}.ubound(1)
                %disp('kww tau over bound');
                newfitinfo.fun{cnt}.curval(1)=newfitinfo.fun{cnt}.ubound(1);
            end
            if newfitinfo.fun{cnt}.curval(1)<newfitinfo.fun{cnt}.lbound(1)
                %disp('kww tau under bound');
                newfitinfo.fun{cnt}.curval(1)=newfitinfo.fun{cnt}.lbound(1);
            end
            if newfitinfo.fun{cnt}.curval(2)>newfitinfo.fun{cnt}.ubound(2)
                %disp('kww beta over bound');
                newfitinfo.fun{cnt}.curval(2)=newfitinfo.fun{cnt}.ubound(2);
            end
            if newfitinfo.fun{cnt}.curval(2)<newfitinfo.fun{cnt}.lbound(2)
                %disp('kww beta under bound');
                newfitinfo.fun{cnt}.curval(2)=newfitinfo.fun{cnt}.lbound(2);
            end        
        end
    end
    % next check t0 bounds
    if newfitinfo.t0curval>newfitinfo.t0ubound
        disp('t0 shift upper bound');
        newfitinfo.t0curval=newfitinfo.t0ubound;
    end
    if newfitinfo.t0curval<newfitinfo.t0lbound
        disp('t0 shift lower bound');
        newfitinfo.t0curval=newfitinfo.t0lbound;
    end
    
    % next check if all parameter values are equal to previous fitinfo, and if so "break"
    % or check if time-limit reached (maybe 5 seconds)
    if cputime-starttime>60
        disp('time-out 60 seconds');
        break;
    end
    %swap in new values for fitinfo
    fitinfo=newfitinfo;
end


% step 2: 10x smaller steps


