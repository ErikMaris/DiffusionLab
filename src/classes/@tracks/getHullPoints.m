function hPs = getHullPoints(XY,n)
% getHullPoints Gets the convex hull points
% 
% -------------------------------------------------------------
% -                         USAGE                             -
% -------------------------------------------------------------
% 
% hPs = getHullPoints(XY,n)

% -------------------------------------------------------------
% -                         HISTORY                           -
% -------------------------------------------------------------
% 
% Reproduced from Florian Meirer 2016
% 
% -------------------------------------------------------------
% Copyright (C) 2021 J.J. Erik Maris
% Inorganic Chemistry & Catalysis, Utrecht University
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

K = convhull(XY(:,1),XY(:,2));
% The convex hull K is expressed in terms of a vector of point indices
% arranged in a counterclockwise cycle around the hull.

%figure(); plot(XY(K,1),XY(K,2),'-x');

% expand the convhull by a value n 
K2 = zeros(size(K));
K2(2:end) = K(1:end-1);
K2(1)=K(end);
v(:,1) = XY(K,1)-XY(K2,1); % each row is the vector between 2 hull points 
v(:,2) = XY(K,2)-XY(K2,2);
% make it length 1:
vl = sqrt(v(:,1).^2 + v(:,2).^2);
v = bsxfun(@rdivide,v,vl);
% remove the NaN from dividing by 0 (lenght of point to itself)
v(isnan(v(:,1)),1)=0;
v(isnan(v(:,2)),2)=0;
% get both normal vectors of each v:
n1(:,1) = -v(:,2);
n1(:,2) = v(:,1);
n2(:,1) = v(:,2);
n2(:,2) = -v(:,1);
% shift normal vectors by 1 point:
n1_shift = zeros(size(n1));
n1_shift(2:end,:)=n1(1:end-1,:);
n1_shift(1,:)=n1(end,:);

n2_shift = zeros(size(n2));
n2_shift(2:end,:)=n2(1:end-1,:);
n2_shift(1,:)=n2(end,:);

% for each point we get 4 new points:
XYnew1(:,1) = XY(K,1) + n1(:,1).*n;
XYnew1(:,2) = XY(K,2) + n1(:,2).*n;

XYnew2(:,1) = XY(K,1) + n2(:,1).*n;
XYnew2(:,2) = XY(K,2) + n2(:,2).*n;

XYnew3(:,1) = XY(K,1) +  n1_shift(:,1).*n;
XYnew3(:,2) = XY(K,2) +  n1_shift(:,2).*n;

XYnew4(:,1) = XY(K,1) + n2_shift(:,1).*n;
XYnew4(:,2) = XY(K,2) + n2_shift(:,2).*n;

XYnew = vertcat(XYnew1,XYnew2,XYnew3,XYnew4);

% check which ones are inside the hull:
in = inpolygon(XYnew(:,1),XYnew(:,2),XY(K,1),XY(K,2));

% new set of points:
XYexthull = XYnew(~in,:);
Knew = convhull(XYexthull(:,1),XYexthull(:,2));
hPs(:,1) = XYexthull(Knew,1);
hPs(:,2) = XYexthull(Knew,2);

% interpolate more points along the hull
dummy=[];
for ii=1:size(hPs,1)-1
    % from first to second point:
    x = [hPs(ii,1) hPs(ii+1,1)]; % two points
    y = [hPs(ii,2) hPs(ii+1,2)];
    if hPs(ii,1)<hPs(ii+1,1)
        xq = (hPs(ii,1):hPs(ii+1,1))'; % x in single pixel steps
        if isempty(xq) || size(xq,1)<10
            xq = (hPs(ii,1):0.1:hPs(ii+1,1))'; % x in finer pixel steps
        end
    else
        xq = (hPs(ii+1,1):hPs(ii,1))';
        if isempty(xq) || size(xq,1)<10
            xq = (hPs(ii+1,1):0.1:hPs(ii,1))'; % x in finer pixel steps
        end
    end    
    yq = interp1(x,y,xq); % interpolate points between the 2 points
    if hPs(ii,1)==hPs(ii+1,1) % x-coordinates are the same
        if hPs(ii+1,2)>hPs(ii,2)
            yq = (hPs(ii,2):hPs(ii+1,2))';
            if isempty(yq) || size(yq,1)<10
                yq = (hPs(ii,2):0.1:hPs(ii+1,2))'; % in finer pixel steps
            end
        else
            yq = (hPs(ii+1,2):hPs(ii,2))';
            if isempty(yq) || size(yq,1)<10
                yq = (hPs(ii+1,2):0.1:hPs(ii,2))'; % in finer pixel steps
            end
        end
        xq = repmat(xq,size(yq,1),1);
    end
    dummy = vertcat(dummy,[xq yq]);
%     hold on; plot(xq,yq,'g-o');
%     drawnow;
%     waitforbuttonpress;
end
clearvars hPs
hPs = unique(dummy,'rows','stable');
% hold on; plot(hPs(:,1),hPs(:,2),'r-o');
% axis image;

% printCurrentFig:
% I moved this to a seperate function m file to be able to call it from
% within script2.m