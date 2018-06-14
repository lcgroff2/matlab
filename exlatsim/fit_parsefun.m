
function fit=fit_parsefun(fit)

% script to minimally parse a function string in order to extract variable names
% and generate a string for fitting purposes

% will be very inefficient, but is intended to be run once, so it's OK

funstr=fit.fnstr;
depvbl=fit.fnxstr;

% first strip out blank spaces (using ismember, but could use isspace)
funstr(ismember(funstr,' '))=[];

% maybe adding one/two blank spaces at the end helps with the logic below
funstr(end+1)=' ';
funstr(end+1)=' ';


% basic rule is to first find an alpha character, then continue until the next non-alphanumeric
% character, and the string between these two points is the variable or function name
curpos=1;
outstr=[];
curvar=0;
parmlist={};
founddependendentvariable=0;

while (curpos+1) < length(funstr),

  % if character at curpos is not alphabetical, move to next character
  if not(isstrprop(funstr(curpos),'alpha')),
    outstr=[outstr funstr(curpos)];
    curpos=curpos+1;
    varnam=''; % this line probably not needed
  else % is alpha, so is the start of a variable name or function name
    % first identify all alphanumeric characters in remainder of string
    isal=isstrprop(funstr((curpos+1):end),'alphanum');
    % then also add in underscore, since those are also valid in variable names
    isal=isal+ismember(funstr((curpos+1):end),'_');
    % now find the first  non-alphanum character--that marks the end of our variable name
    nonalphidx=find(not(isal));
    if length(nonalphidx)>0,
      firstnonalph=nonalphidx(1);      
      % extract variable name, and increment curpos to start of next possible variable occurrence
      varnam=funstr(curpos:(curpos+firstnonalph-1));
      curpos=curpos+firstnonalph;
    else  % otherwise is a 1-character variable name
      varnam=funstr(curpos);
      curpos=curpos+1;
    end
    % check to see if it is a function or other keyword
    if funstr(curpos)=='(',
      %fprintf(1,'%s is a function\n',varnam);
      outstr=[outstr varnam];
    elseif strcmp(varnam,depvbl), % check to see if is the dependent variable
      %fprintf(1,'%s is a dependent variable\n',varnam);
      outstr=[outstr fit.xstr]; % insert dependent variable name
      %
      founddependendentvariable=1;
    else
      %fprintf(1,'%s is a parameter\n',varnam);
      % if not, then check against current parameter list.  If not,
	  % add variable name to list, increment parameter number, and add 
	  % "parm(XX)" to outstr, otherwise it's a bit complicated
	  if length(parmlist)==0, % no parameters to check, so add to parmlist
        curvar=curvar+1;
        outstr=[outstr sprintf('parm(%i)',curvar)];
        parmlist{curvar}=varnam;
	  else
	    foundmatch=0;
	    for idx=1:length(parmlist), % look for a match
		  if strcmp(parmlist{idx},varnam), % found a match, so use the same parm#
		    outstr=[outstr sprintf('parm(%i)',idx)];
			foundmatch=1;
		  end
		end
		if not(foundmatch), % no match, so add to parmlist, outstr as normal
		  curvar=curvar+1;
		  outstr=[outstr sprintf('parm(%i)',curvar)];
		  parmlist{curvar}=varnam;
		end
	  end
    end
  end % end if
  %outstr;
  
end % end while



if not(founddependendentvariable),
  disp('warning: no dependent variable found')
end

fit.parsedfnstr=outstr;
fit.parmlist=parmlist;
