% Function compute NC_tune, NB_tune
function [lambda_x] = len_param(X,lambda,M_Y)
if ~isempty(X)
    Mat_X = [];
    for i = 1:size(X,1)
        tmp =[];
        for j = i+1:size(X,1)
            tmp = [tmp (norm(X(i,:)-X(j,:)))^2];
        end
        Mat_X = [Mat_X tmp];
    end
    M_X = median(Mat_X);
    lambda_x = (lambda*M_Y)/(M_X);
else
    lambda_x = [];
end
end
