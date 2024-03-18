function error_returned = errorMap(target, ML, errorType)
orishape = cat(3, ones(16,1)*linspace(0,40,16), linspace(0,40,16)'*ones(1,16), zeros(16,16));
disp = sqrt(sum((target-orishape).^2,3));
dist = sqrt(sum((target-ML).^2,3));

if strcmp(errorType,'absolute')
    error_returned = dist;

    figure(); 
    s = pcolor(error_returned); 
    cb=colorbar; 
    cb.Label.String='Error(mm)';
else
    if strcmp(errorType,'relative1')
        error_returned = dist / max( disp(:) );
    elseif strcmp(errorType,'relative2')
        error_returned = dist / 40;
    end

    figure(); 
    s = pcolor(error_returned*100); 
    cb=colorbar; 
    cb.Label.String='Error(%)';
end
s.FaceColor = 'interp';
s.EdgeColor = 'w';
title(['Max ',num2str(max(error_returned(:)*100),3),'%, Avg ',num2str(mean(error_returned(:)*100),3),'%']);
set(gca,'FontName','Arial','FontSize',15);
colormap(viridis);
end