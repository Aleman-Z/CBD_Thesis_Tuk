function out = normalize_matrix(mat, N)

rem_norm = reshape(mat, 1, []);
rem_norm = rem_norm(reshape(N,1,[]));
rem_norm = normalize(rem_norm);
mat2 = mat;
mat2 = reshape(mat2,1,[]);
mat2(reshape(N,1,[])==1) = rem_norm;
mat2 = reshape(mat2, size(mat,1),size(mat,2));
out = mat2;
end

