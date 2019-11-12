function CH = computeCH(Aeq,lbs,ubs,dims,np,nc,tolz,tolp)
    % Computes the convex hull for production envelopes of metabolic network. Solution is the list of hyperplanes and
    % set of extreme points of the Convex hull. Inputs are:
    % Aeq: the Stoichiometric matrix
    % lbs ubs : and constraints on the fluxes (upper and lower bounds)
    % dims: indices of dimensions onto which the projection should be computed
    % Optional arguments:
    % n_dec_p: general precision for arithmetic operations
    % (number of decimals-integer value). Default is 12.
    % n_dec_c: precision for value comparison >=,<=,=,etc. Default is 6.
    % tol_zero: value below which a variable is considered to be zero
    % tolp: when adding an equality constraint (Ax = b), if infeasible,
    % add Ax <= b + tolp ,  Ax >= b - tolp   Default is 1e-5

    % computeCH returns 2 variables:
    % - CH.eps contains the extreme points
    % - CH.hps contains the hyperplanes

    global chull;
    global ePoints;
    global tol_zero;
    global tol_mem;
    global n_dec_c;
    global n_dec_p;
    global tol_lp;

    if(nargin<4)
        error('Error using computeCH. Not enough input arguments.')
    end
    if(nargin == 4)
	% default values
        n_dec_p = 12;
        n_dec_c = 6;
        tol_mem = 1*10^-n_dec_c; % format used by the function uniquetol of matlab
        tol_zero = 1e-5;
        tol_lp = 1e-5;
    else
        n_dec_c = nc;
        n_dec_p = np;
        tol_mem = 1*10^-n_dec_c;
        tol_zero = tolz;
        tol_lp = tolp;
    end

    num_mets= size(Aeq,1);
    beq = zeros(num_mets,1);
    A = [];
    b = [];

    polytope = {A,b,Aeq,beq,lbs,ubs,dims};
    % Initial points
    ePoints  = InitialPoints(polytope);
    % Initial Hull
    chull = InitialHull(ePoints,dims);
    % Hull refinement
    IncrementalRefinement(polytope);
    % return list of HPs
    ht = GetListHPs(dims);
    CH ={};
    CH.hps = ht;
    CH.eps = ePoints(dims,:).';
end
