// Velocity for regularization term.
// Computed by filtering the exactly deconvoluted velocity field V
// with filterWidthRatioReg times the base filter width.
solve(fvm::ddt(W) + fvm::Sp(1/filterWidthReg,W) == 1/filterWidthReg*V);



