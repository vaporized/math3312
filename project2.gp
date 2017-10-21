/* GP Code for MATH3312 Project 2
 * author: vapor
 */

 \\convert fraction to real
n(x)=x/1.;

divided_diff_pol(xs, ys)={
	my(F=Map(), len, coefs, tmp, ret=0);
	\\Note that len == n+1 in the discription of algorithm
	if((len=length(xs))!=length(ys), return([]););
	for(i=0,len-1,mapput(F,[i,0],ys[i+1]));
	for(i=1,len-1,
		for(j=1,i,
			mapput(F,[i,j],(mapget(F,[i,j-1])-mapget(F,[i-1,j-1]))/(xs[i+1]-xs[i-j+1]));
		);
	);
	for(i=0,len-1,
		tmp=1;
		for(j=1,i,tmp*=x-xs[j];);
		ret+=mapget(F,[i,i])*tmp;
	);
	return(ret);
}

pol_eval(pol, val)={
	my(f(x)=eval(pol));
	return(n(f(val)));
}


\\main
years=[1940,1950,1960,1970,1980,1990];
populations=[132165,151326,179323,203302,226542,249633];

default(format,"g0.12");

est_pol=divided_diff_pol(years,populations);
print("The population in the year 1930 is: ",pol_eval(est_pol,1930));
print("The population in the year 1965 is: ",pol_eval(est_pol,1965));
print("The population in the year 2000 is: ",pol_eval(est_pol,2000));
print("The fitted polynomial is: ", n(est_pol));

quit;



