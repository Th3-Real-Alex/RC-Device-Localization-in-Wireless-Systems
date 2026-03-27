function helperPlotTOADevicePositions(c,toaEst,tgtPosEst,anchorPos,tgtPos)
%helperPlotTOADevicePositions Plots 2D positions of all anchors and
%targets, estimated positions of targets, and trilateration circles

% Range estimation
rngEst = toaEst*c;

% Plot positions of all anchors, target, and estimated position of the target
figure
plot(anchorPos(1,:),anchorPos(2,:),'b^','LineWidth',2,'MarkerSize',10),hold on
plot(tgtPos(1,:),tgtPos(2,:),'rx','LineWidth',2,'MarkerSize',10),hold on
plot(tgtPosEst(1),tgtPosEst(2),'go','LineWidth',2,'MarkerSize',12),hold on

% Get distance from anchors to target to keep consistent plot axes
radiustx = sqrt(anchorPos(1,:).^2 + anchorPos(2,:).^2);
range = max(abs(anchorPos+repmat(radiustx,[size(anchorPos,1) 1])),[],'all');
xlim([-range range]),ylim([-range range]);
axis equal
grid on;
xlabel('x-axis (meters)'),ylabel('y-axis (meters)');

% Plot trilateration circles
numAnchor = length(rngEst);
angles = 0:2*pi/720:2*pi;
for anchorIdx = 1:numAnchor
    x = rngEst(anchorIdx) * cos(angles) + anchorPos(1,anchorIdx);
    y = rngEst(anchorIdx) * sin(angles) + anchorPos(2,anchorIdx);
    plot(x,y,'c--','LineWidth',1);
    hold on; grid on
end
legend({'Anchor Positions','Target Position','TOA Position Estimate','Trilateration Circles'},'Location','best','FontSize',10)
end