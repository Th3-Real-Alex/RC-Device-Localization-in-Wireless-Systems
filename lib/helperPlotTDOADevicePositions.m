function helperPlotTDOADevicePositions(c,tdoaEst,tgtPosEst,anchorPos,tgtPos,varargin)
%helperPlotTDOADevicePositions Plots 2D positions of all anchors and
%targets, estimated positions of targets, and hyperbola curves

% Range difference estimation
rngDiffEst = tdoaEst*c;

% Plot positions of all anchors, target, and estimated position of the target
figure
plot(anchorPos(1,:),anchorPos(2,:),'b^','LineWidth',2,'MarkerSize',10), hold on
plot(tgtPos(1,:),tgtPos(2,:),'rx','LineWidth',2,'MarkerSize',10), hold on
if nargin>5
    txPos = varargin{1};
    plot(txPos(1,:),txPos(2,:),'kv','LineWidth',2,'MarkerSize',10), hold on
end
plot(tgtPosEst(1),tgtPosEst(2),'go','LineWidth',2,'MarkerSize',10), hold on

% Plot hyperbola curves
numAnchorpair = length(rngDiffEst);
for anchorPairIdx = 1:numAnchorpair
    % Get hyperbolic curve for the TDOA between an anchor
    % and anchor 1
    [x, y] = get2DHyperbolicSurface( ...
        anchorPos(:,1), ...
        anchorPos(:,anchorPairIdx+1), ...
        rngDiffEst(anchorPairIdx));
    if isreal(x) && isreal(y)
        plot(x,y,'c--','LineWidth',1)
    end
    hold on
end

% Get distance from APs to STA to keep consistent plot axes
radiustx = sqrt(anchorPos(1,:).^2 + anchorPos(2,:).^2);
range = max(abs(anchorPos+repmat(radiustx,[size(anchorPos,1) 1])),[],'all');
xlim([-range range]),ylim([-range range]);
grid on;
xlabel('x-axis (meters)'),ylabel('y-axis (meters)');
if nargin>5
    legend({'Anchor Positions','Target Position','Separate Transmitter Position','TDOA Position Estimate','Hyperbola Curves'}, ...
        'Location','best','FontSize',10);
else
    legend({'Anchor Positions','Target Position','TDOA Position Estimate','Hyperbola Curves'}, ...
        'Location','best','FontSize',10);
end
end

function [x, y] = get2DHyperbolicSurface(anchorRefPos,anchorPos,rngDiffEst)
%get2DHyperbolicSurface Get 2D hyperbolic surface for a given pair of anchors
%
%   [X,Y] = get2DHyperbolicSurface(ANCHORREFPOS,ANCHORPOS,RNGDIFFEST) finds
%   a 2D hyperbolic surface between two anchors, in a manner that all
%   points of the surface correspond to the same TDOA. The input
%   ANCHORREFPOS is the reference anchor postion, ANCHORPOS is another
%   anchor position, and RNGDIFFEST is the Difference of target ranges
%   between the anchor pair and the target. The output X and Y are x-y
%   coordinates of the hyperbolic surface.

theta = linspace(-pi/2,pi/2,100);
phi = 0;
[Theta, Phi] = meshgrid(theta,phi);

% Compute hyperbola quantities
D = norm(anchorRefPos - anchorPos)/2;
c = rngDiffEst;

% Hyperbola equation in a frame where sensors lie on x axis and are located
% at -D and D on x axis
xC = -c./cos(Theta)./2;
yC = sqrt(4*D^2 - c^2).*tan(Theta).*cos(Phi)./2;
zC = sqrt(4*D^2 - c^2).*tan(Theta).*sin(Phi)./2;

% Translate and rotate to the scenario frame
r0 = (anchorPos + anchorRefPos)/2;
a = [1;0;0];
b = (anchorPos - anchorRefPos);
b = b/norm(b);
v = cross(a,b);
s = norm(v);
c = dot(a,b);
V = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
if abs(s) > 0
    R = eye(3) + V + V^2*(1 - c)/s^2;
else
    R = eye(3);
end

x = R(1,1).*xC + R(1,2).*yC + R(1,3).*zC;
y = R(2,1).*xC + R(2,2).*yC + R(2,3).*zC;
x = x + r0(1);
y = y + r0(2);
end