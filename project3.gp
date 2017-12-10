/* GP Code for MATH3312 Project 3
 * author: vapor
 */

\\convert fraction to real
n(x)=x*1.;

 \\get diagonal component of a square matrix and its remainder
 matdiag(m) = {
 	my(identity_mat, diag, rem);
 	identity_mat = matid(#m);
 	diag = matrix(#m[,1],#m,i,j,m[i,j]*identity_mat[i,j]);
 	rem = m - diag;
 	return([diag,rem]);
 }

jacobi_iter(A, b, init_guess, max_iter = 30, tol = 0.001) = {
	my(d, r, res = init_guess, prev);
	[d, r] = matdiag(A);
	for(iter=1,max_iter,
		prev = res; 
		res = 1/d * (b - r*prev);
		if(vecmax(abs(res - prev))<=tol, 
			return(res);
		);
	);
	print("Failed to converge within ", max_iter, " loops.");
	return(res);
}

input_A = [10, -1, 2, 0; -1, 11, -1, 3; 2, -1, 10, -1; 0, 3, -1, 8];
input_b = [6, 25, -11, 15]~;
input_guess = [0, 0, 0, 0]~;

print("Input Matrix A: ", input_A);
print("Input Vector b: ", input_b);
print("Answer from Jacobi Iterative Method: ", n(jacobi_iter(input_A, input_b, input_guess)));
print("Exact solution: ", matsolve(input_A, input_b));

quit;