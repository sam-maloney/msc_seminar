
// adjust temporal stress tensor tau for the next time step
	
solve(fvm::ddt(tau) + fvm::Sp(1/filterWidth,tau) == filterWidth*fvc::ddt(U)*(fvc::ddt(U)));		

// calculate V

V = U + filterWidth * fvc::ddt(U);
