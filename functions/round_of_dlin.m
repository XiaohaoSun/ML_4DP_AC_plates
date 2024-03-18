function dlin_round = round_of_dlin(dlin)
% falls in range (-1,1)
dlin_opti = tanh(dlin);
dlin_round = dlin_opti;
dlin_round(dlin_round>0) = 1; dlin_round(dlin_round<0) = -1; 
dlin_round = gather(extractdata(dlin_round));
end
