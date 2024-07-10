function x = rejection_sampling(fun,maxValue)
accepted = 0;
while accepted == 0
    % Proposal value
    proposal = rand(1);
    % Function value on propsal point
    funValue = fun(proposal);
    vertValue = rand(1)*maxValue;
    if vertValue < funValue
        x = proposal;
        accepted = 1;
    end
end
end