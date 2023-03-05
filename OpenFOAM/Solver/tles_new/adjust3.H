// clean divergence of U using projection scheme

solve(fvm::laplacian(U_) == fvc::div(U));
U = U - fvc::grad(U_);
