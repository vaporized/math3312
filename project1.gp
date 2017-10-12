/* GP Code for MATH3312 Project 1
 * author: vapor
 *
 * The script includes both iterative and recursive implementations for the two methods.
 * The iterative ones are from the textbook, the recursive ones are more concise and have no side effeccts.
 * The answer to task 3 is obvious: bisection converges in 32 iterations and the mixed one in 9, which indicates
 *  Newton's method is better for this example.
 */

\\convert fraction to real
n(x)=x+0.;


bisection_iterative(fun,lower,upper,tol=10E-10,maxiter=10E3,verbose=0)={
	my(mid,flow=fun(lower),fmid,iter=1);
	while(iter<=maxiter,
		mid=lower+(upper-lower)/2;
		fmid=fun(mid);
		if(verbose,print("p",iter," is: ",mid+0.));
		if(fmid==0||(upper-lower)/2<tol,
			return([mid,iter]);
		);
		iter++;
		if(flow*fmid>0,
			lower=mid;flow=fmid,
			upper=mid);
		);
		print("Failed to converge within ",maxiter," iterations");
	return([mid,iter-1])
}

\\index starts at 1 because the initial average is p1
bisection_recursive(fun,lower,upper,tol=10E-10,maxiter=10E3,curiter=1,verbose=0)={
	my(mid=lower+(upper-lower)/2,flow=fun(lower),fmid=fun(mid));
	if(verbose,print("p",curiter," is: ",mid+0.));
	if(maxiter==1,
		print("Failed to converge within ",curiter," iterations");
		return([mid,curiter]);
	);
	if(fmid==0||(upper-lower)/2<tol,
		return([mid,curiter]);
	);
	if(flow*fmid>0,
		return(bisection_recursive(fun,mid,upper,tol,maxiter-1,curiter+1,verbose)),
		return(bisection_recursive(fun,lower,mid,tol,maxiter-1,curiter+1,verbose)))
}

newton_iterative(fun,guess,tol=10E-10,maxiter=10E3,verbose=0)={
	my(iter=1,p,df=fun');
	while(iter<=maxiter,
		p=guess-fun(guess)/df(guess);
		if(verbose,print("p",iter," is: ",p+0.));
		if(abs(p-guess)<=tol,return([p,iter]););
		iter++;guess=p;
	);
	print("Failed to converge within ",maxiter," iterations");
	return([guess,iter-1])
}

\\index starts at 0 because the initial guess is p0
newton_recursive(fun,guess,tol=10E-10,maxiter=10E3,curiter=0,verbose=0)={
	my(p=guess-fun(guess)/fun'(guess));
	if(maxiter==0,print("Failed to converge within ",curiter," iterations");return([guess,curiter]));
	if(abs(p-guess)<=tol,
		return([p,curiter+1]),
		return(newton_recursive(fun,p,tol,maxiter-1,curiter+1,verbose)))
}

f(x)=-x^3-cos(x);

task1_i=bisection_iterative(f,-2,2);
task1_r=bisection_recursive(f,-2,2);

print("task1(iterative): ",n(task1_i[1]),", ",task1_i[2],"iterations.");
print("task1(recursive): ",n(task1_r[1]),", ",task1_r[2],"iterations.");

task21_i=bisection_iterative(f,-2,2,10E-10,4);
task21_r=bisection_recursive(f,-2,2,10E-10,4,1);

print("task2-part1(iterative): ",n(task21_i[1]),", ",task21_i[2],"iterations.");
print("task2-part1(recursive): ",n(task21_r[1]),", ",task21_r[2],"iterations.");

newton_init=task21_i[1];
print("value: ",n(newton_init)," is used for the initial guess");

task22_i=newton_iterative(f,newton_init);
task22_r=newton_recursive(f,newton_init);

print("task2-part2(iterative): ",n(task22_i[1]),", ",task22_i[2],"iterations.");
print("task2-part2(recursive): ",n(task22_r[1]),", ",task22_r[2],"iterations.");

\\task3_i=newton_iterative(f,newton_init,10E-10,4);
\\task3_r=newton_recursive(f,newton_init,10E-10,4,0);

\\print("task3(iterative): ",n(task3_i[1]),", ",task3_i[2],"iterations.");
\\print("task3(recursive): ",n(task3_r[1]),", ",task3_r[2],"iterations.");