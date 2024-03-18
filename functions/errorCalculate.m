function error_returned = errorCalculate(target, ML, errorType)
orishape = cat(3, ones(16,1)*linspace(0,40,16), linspace(0,40,16)'*ones(1,16), zeros(16,16));
disp = sqrt(sum((target-orishape).^2,3));
dist = sqrt(sum((target-ML).^2,3));

if strcmp(errorType,'absolute')
    error_returned = dist;
elseif strcmp(errorType,'relative1')
    error_returned = dist / max( disp(:) );
elseif strcmp(errorType,'relative2')
    error_returned = dist / 40;
end

end