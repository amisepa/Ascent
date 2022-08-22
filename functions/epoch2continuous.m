% Converts an epoched dataset into a continuous one.
%          Data segments are concatenated using a 'boundary' events.
%
% Usage: EEG = epoch2continuous(EEG);
%
% Cedric Cannard, August 2022

function EEG = epoch2continuous(EEG)

if nargin<1, help epoch2continuous; return; end
if isempty(EEG.epoch), error('epoch2continuous() only works for epoched data!'); end

[~,indx,~] = unique([EEG.event.latency],'first','legacy');
nEvents  = length(indx);

% new type
typearray  = {EEG.event(indx).type};

% new duration
if isfield(EEG.event, 'duration')
      durarray = {EEG.event(indx).duration};
else
      durarray = num2cell(ones(1,nEvents));
end
if ~isfield(EEG.event, 'urevent') 
      [EEG.event.urevent] = EEG.event.type;
end

% new urevent
evArray = {EEG.event(indx).urevent};
nEpochs = EEG.trials;
nPnts = EEG.pnts;
latsamarray = zeros(1, nEvents);

% new continuous latencies
for iEv=1:nEvents
    ep  = EEG.event(indx(iEv)).epoch;
    if iscell(EEG.epoch(ep).eventlatency)
        target = EEG.epoch(ep).eventlatency{EEG.epoch(ep).event == indx(iEv)};
    else
        target = EEG.epoch(ep).eventlatency;
    end
    times = EEG.times;
    narginchk(1,2)
    nargoutchk(0,3)
    ntarget = length(target);
    [cvalue, cindex, cdiff] = deal(zeros(1, ntarget));
    for iTarget = 1:ntarget
        if isnan(target(iTarget)) || isinf(target(iTarget))
            [cvalue(iTarget), cindex(iTarget), cdiff(iTarget)] = deal(NaN);
        else
            [~, cindex(iTarget)] = min(abs(times-target(iTarget)));
            cvalue(iTarget) = times(cindex(iTarget));
            cdiff(iTarget)  = cvalue(iTarget)-target(iTarget); % times minus lat
        end
    end
    latsamarray(iEv) = cindex + (ep-1)*nPnts;
end
latsamarray = num2cell(latsamarray);

% new boundaries
latbound    = num2cell(nPnts:nPnts:nPnts*nEpochs);
nbound      = length(latbound);
boundarray  = repmat({'boundary'}, 1, nbound);

% concatenates events and boundaries info
typearray   = [typearray boundarray];
latsamarray = [latsamarray latbound];
evArray     = [evArray repmat({0},1,nbound)];
durarray    = [durarray repmat({0},1,nbound)];
nEvents   = length(typearray); % new length

% Builts new EEG
EEG.trials  = 1;
EEG.xmin    = 0;
EEG.data    = reshape(EEG.data , EEG.nbchan,nEpochs*nPnts);
EEG.event   = [];

% Events
[EEG.event(1:nEvents).type    ] = typearray{:};
[EEG.event(1:nEvents).latency ] = latsamarray{:};
[EEG.event(1:nEvents).urevent ] = evArray{:};
[EEG.event(1:nEvents).duration] = durarray{:};
EEG.epoch  = [];
EEG.times  = [];
EEG.epochdescription = {};
EEG.reject = [];
EEG.pnts   = length(EEG.data);
EEG.xmax   = length(EEG.data)/EEG.srate;

EEG = eeg_checkset(EEG, 'eventconsistency');
