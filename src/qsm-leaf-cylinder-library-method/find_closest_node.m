function [V,I] = find_closest_node(nodes,value)

% Find the node closest to the given value

if isscalar(nodes)
    I = 1;
elseif value > max(nodes)
    I = length(nodes);
else
    % First node above or equal to value
    iNodeAbove = find(nodes >= value,1,'first');
    if iNodeAbove == 1
        I = iNodeAbove;
    else
        % First node below value
        iNodeBelow = iNodeAbove - 1;
        % Choose the node closer to value
        if abs(value-nodes(iNodeBelow)) < abs(value-nodes(iNodeAbove))
            I = iNodeBelow;
        else
            I = iNodeAbove;
        end
    end
end
V = nodes(I);

end