
function beta = linfit(x,y)

beta0(:,1:10) = 1;

beta = nlinfit(x,y,@model,beta0);
