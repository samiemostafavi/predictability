function figCoord = axpos2figpos(ax, dataCoord)
    % Convert axes data point to normalized figure position
    axUnits = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    axPos = get(ax, 'Position');
    axLimX = xlim(ax);
    axLimY = ylim(ax);
    
    xRel = (dataCoord(1) - axLimX(1)) / diff(axLimX);
    yRel = (dataCoord(2) - axLimY(1)) / diff(axLimY);
    
    figCoord = [axPos(1) + axPos(3) * xRel, ...
                axPos(2) + axPos(4) * yRel];
    set(ax, 'Units', axUnits);  % restore original
end