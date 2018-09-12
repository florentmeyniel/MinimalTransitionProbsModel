function [NBgA, NAgB, NAgA, NBgB, s1, NA, NB] = CountEventInMemory(s, param)
% Compute transitions in sequence s.
%
% Usage: [NBgA, NAgB, NAgA, NBgB, s1, NA, NB] = CountEventInMemory(s, param)
%   Input:
%           s: sequence of event (1s and 2s)
%           param: cell with paired arguments (name, value)
%               * the default is to use all events without decay
%               * {'Limited', 10}: the memory window size is limited to 10 events
%               * {'Decay', 2} an exponential decay exp(-n/2) is applied 
%                 to event n in the past. Note: the latest even is n=1
%               * can be combined, e.g. {'Limited', 10, 'Decay', 2}
%
% Output:
%           Nigj: transition count from j to i
%           Ni: event count for event type i
%           s1: 1st event of the sequence considered.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
windowed = 0; MemDecay = Inf;
if exist('param', 'var')
    if iscell(param) && any(length(param) == [2 4])
        for k = find(mod(1:length(param), 2)) % take only odd indices
            if strcmp(param{k}, 'Limited')
                windowed = 1;
                MemSpan = param{k+1};
            end
            
            if strcmp(param{k}, 'Decay')
                MemDecay = param{k+1};
            end
        end
    elseif ~isempty(param)
        error('check parameters')
    end
end

% Compute the windowed sequence to consider
if windowed == 1 && length(s) >= MemSpan
    % get only recent events
    trn      = diff(s(end-MemSpan+1:end)); % 1: A -> B; -1: B -> A;
    subs4trn = s(end-MemSpan+1:end-1);
    subs     = s(end-MemSpan+1:end);
    
    % 1st event
    s1 = s(end-MemSpan+1);
else
    % get all event if they are fewer events than the memory span
    trn      = diff(s(1:end)); % 1: A -> B; -1: B -> A;
    subs4trn = s(1:end-1);
    subs     = s;
    
    % 1st event
    s1 = s(1);
end

if MemDecay < Inf
    if nargout <=5
        % the caller requests count about transitions
        % for speed, compute only that
        
        % Compute the decay factor
        TrnDecay = exp(-(1/MemDecay)*(length(subs4trn)+1 - [1:length(subs4trn)]) );
        
        % compute observed transition (with decay & memory span)
        NBgA = sum( (trn(subs4trn==1)==1)  .* TrnDecay(subs4trn==1));
        NAgB = sum( (trn(subs4trn==2)==-1) .* TrnDecay(subs4trn==2));
        NAgA = sum( (trn(subs4trn==1)==0)  .* TrnDecay(subs4trn==1));
        NBgB = sum( (trn(subs4trn==2)==0)  .* TrnDecay(subs4trn==2));
        
    elseif nargout == 2
        % the caller requests count about stimuli
        % for speed, compute only that
        
        % Compute the decay factor
        SDecay   = exp(-(1/MemDecay)*(length(subs)    +1 - [1:length(subs)    ]) );
        
        % compute observed event count (with decay & memory span)
        NA = sum(SDecay(subs==1));
        NB = sum(SDecay(subs==2));
        
    else
        % the call request everything
        
        % Compute the decay factor
        TrnDecay = exp(-(1/MemDecay)*(length(subs4trn)+1 - [1:length(subs4trn)]) );
        SDecay   = exp(-(1/MemDecay)*(length(subs)    +1 - [1:length(subs)    ]) );
        
        % compute observed transition (with decay & memory span)
        NBgA = sum( (trn(subs4trn==1)==1)  .* TrnDecay(subs4trn==1));
        NAgB = sum( (trn(subs4trn==2)==-1) .* TrnDecay(subs4trn==2));
        NAgA = sum( (trn(subs4trn==1)==0)  .* TrnDecay(subs4trn==1));
        NBgB = sum( (trn(subs4trn==2)==0)  .* TrnDecay(subs4trn==2));
        
        % compute observed event count (with decay & memory span)
        NA = sum(SDecay(subs==1));
        NB = sum(SDecay(subs==2));
    end
    
else
    % although the previous lines are correct whatever MemDecay, for speed
    % avoid computing decay when there is not (= always 1)
    % compute observed transitions (within memory span)
    NBgA = sum( trn(subs4trn==1)==1 );
    NAgB = sum( trn(subs4trn==2)==-1);
    NAgA = sum( trn(subs4trn==1)==0 );
    NBgB = sum( trn(subs4trn==2)==0 );
    
    % compute only if asked
    if nargout == 7 || nargout == 2
        NA = sum(subs==1);
        NB = sum(subs==2);
    end
end

end
