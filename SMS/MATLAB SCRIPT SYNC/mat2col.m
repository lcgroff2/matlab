function [v,dim]=mat2col(z)

% MAT2COL -- flattens a 2D matrix to a column vector.  Useful for
% functions that require a vector, such as using mldivide for
% linear least squares involving matrix data (e.g., fitting a 3D
% gaussian).

dim=size(z);
vlen=dim(1)*dim(2);
v=reshape(z,[vlen 1]);