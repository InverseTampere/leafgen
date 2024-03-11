function [V,I] = find_closest_node(nodes,value)

% Ensure the value is between minimum and maximum node
if value < min(nodes) || value > max(nodes)
    error('Input value is outside the node range.')
end    

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
V = nodes(I);
end