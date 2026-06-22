%% function
function photons=brightness(Ion,Ton,Toff,Tbl,frames)
%Simulate the intensity trace of an emitter (photons per frame).
%
%Inputs:
% Ion       maximum signal per emitter per frame [photons]
% Ton       average duration of the on-state [frames]
% Toff      average duration of the off-state [frames]
% Tbl       bleaching lifetime [frames]
% frames    number of frames comprising the image sequence [frames]
%
%Outputs:
% photons   intensity trace of an emitter [photons]

cycle=Ton + Toff; % length of a cycle: for a fluorophore to reach the onstate and in the offstate 
cycles=10 + ceil(frames/cycle); % number of cycles in the entire experiment (the +10 creates ten cycles in addition to avoid problems near t~0)
times=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))]; % probability of being in one state is between 0 and 1 (normal distribution)
times(1)=times(1) - rand*(sum(times(1:10))); % for t~0, the probability is one. Therefore, this makes the fluorophore start in the state Toff (in theory)
% otherwise the initial values are far to high
times=cumsum(times(:)); % cumulative sums of times
% --- redo the exact same steps if the times is not as long as frames
while times(end) < frames
   cycles=ceil(2*(frames - times(end))/cycle);
   cycles=[-Toff*log(rand(1,cycles));-Ton*log(rand(1,cycles))];
   cycles(1)=cycles(1) + times(end);
%    temp = elongate(times,cycles(:));clear times;
%    times = temp; clear temp;
   times=[times;cumsum(cycles(:))];
end

times=times.';
% times contains successively periods of activation and periods of off
% (hence the blinking) both described by a normal distribution. Since Toff
% > Ton, then the fluorophore is longer in the off state.

Ton=times(2:2:end) - times(1:2:end); % times(2,:) when size(times) was 2x60 (line 46)
Tbl=cumsum(Ton) + Tbl*log(rand); % the bleaching state is another state as 
% Ton and Toff. Here we cumulate only the times of Ton since it is the one
% that contributes to approach the bleaching state. Tbl follows the same
% distribution as the other two states (hence the log).
n=find(Tbl > 0); % did Ton reached Tbl ? 
if any(n)
   Ton(n(2:end))=0; % the fluorophore is bleached (no Ton left in the signal)
   n=n(1);
   Ton(n)=Ton(n) - Tbl(n);
   times(2*n)=times(2*n) - Tbl(n); % the last "times" is put to 0
end
photons=[zeros(size(Ton));Ion*Ton];clear Ion; clear Ton; clear Tbl;
photons=cumsum(photons(:)); % one point over two has no photons, this is 
% the period in the off state. All others have photons for successive
% duratoins of Ton.
photons=diff(interp1(times,photons,0:frames,'linear',0));clear times;
% times' unit is in number of frames. So here, we are actually looking for
% the number of photons in specific frames (which can in principle not be
% sampled in "times" thus the interpolation). We are computing photons(times)
% at times = 0:frames. Toff > Ton --> thus big steps in times (or large
% periods of frames) corresponds to odd samples (sample point that are
% defined by Toff), and small steps in times (or small periods of frames)
% corresponds to even samples (sample point that are defined by Ton).
% Therefore when looking for the number of photons in a arbitrary frame, it
% is likely the number of fluorophore is small (since there are much more
% frames "stay off" then frames "go on").
end
