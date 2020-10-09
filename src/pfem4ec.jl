"""
pfem4ec

pfem4ec is a fast in-core solver to compute effective electrical conductivity of heterogeneous 
materials from raw images using Pixel-Based Finite Element Method (PFEM). 

See https://github.com/LCC-UFF/pfem4ec/tree/master/docs for documentation.
"""

# Packages
using JSON
using SparseArrays
using LinearAlgebra

# Immutable Model data struct
"""
    Immutable Model data struct
   
    The struct covering the fundamental data from JSON
"""
struct Model
    m_nx::UInt64;
    m_ny::UInt64;
    m_refinement::UInt64;
    m_nMat::UInt16;
    m_rhsType::UInt8;
    m_solverType::UInt8;
    m_pcgTol::Float64;
    m_pcgIter::UInt64;
    m_matKeys::Array{UInt16};
    m_matProp::Array{Float64}; 
    m_nDOFNode::UInt8;
    m_nNodes::UInt64;
    m_nElems::UInt64;
    m_nDOFs::UInt64;
    function Model(_nx,_ny,_refine,_nMat,_rhsType,_solverType,_pcgTol,_pcgIter,_matKeys,_matProp,_nDOFNode,_nNodes,_nElems,_nDOFs)
        new(_nx,_ny,_refine,_nMat,_rhsType,_solverType,_pcgTol,_pcgIter,_matKeys,_matProp,_nDOFNode,_nNodes,_nElems,_nDOFs);
    end
end

# Read JSON file:
"""
    readJSON!(_filename::String)

    Containing specific parameters to guide the FEM analysis. 

    image_dimensions, refinement, type_of_rhs, type_of_solver, solver_tolerance, 
    number_of_iterations, number_of_materials, properties_of_materials
"""
function readJSON!(_filename::String)
    println(".Read JSON!")
    # Open and read file:
    open(_filename,"r") do f
        data = JSON.parse(f);
        nx = data["image_dimensions"][1];
        ny = data["image_dimensions"][2];
        refinement = 1;
        if haskey(data, "refinement"); refinement = data["refinement"]; end
        rhsType = 0;
        if haskey(data, "type_of_rhs"); rhsType = data["type_of_rhs"]; end
        solverType = 0;
        if haskey(data, "type_of_solver"); solverType = data["type_of_solver"]; end
        pcgTol = 0.000001;
        if haskey(data, "solver_tolerance"); pcgTol = data["solver_tolerance"]; end
        pcgIter = nx*ny*refinement*refinement;
        if haskey(data, "number_of_iterations"); pcgIter = data["number_of_iterations"]; end
        nMat = data["number_of_materials"];
        materials = data["properties_of_materials"];
        matKeys = zeros(UInt16,256);
        matProp = zeros(Float64,256);
        for i=1:nMat
            matKeys[convert(UInt8,materials[i][1])+1] = i;
            matProp[convert(UInt8,materials[i][1])+1] = convert(Float64, materials[i][2]);
        end
        # Update the Model based on the given refinement level:
        nDOFNode = 1;
        nx *= refinement;
        ny *= refinement;
        nNodes = (nx+1)*(ny+1);
        nElems = (nx)*(ny);
        nDOFs = nDOFNode*(nx)*(ny);
        # Build the Model with the data read from JSON file:
        model = Model(nx,ny,refinement,nMat,rhsType,solverType,pcgTol,pcgIter,matKeys,matProp,nDOFNode,nNodes,nElems,nDOFs)
        materials = nothing;
        matKeys = nothing;
        matProp = nothing;
        data = nothing;
        return model
    end
end

# Read RAW file:
"""
    readRAW!(_model::Model, _elemMatMap::Array{UInt16,1}, _filename::String)

    The first input is a binary RAW file, that represents images of the microscale 
    of a heterogeneous material with 8-bit data (integer values within the range 0 to 255).     

"""
function readRAW!(_model::Model, _elemMatMap::Array{UInt16,1}, _filename::String)
    println(".Read RAW!");
    # Initializations:
    refine = _model.m_refinement; nx = _model.m_nx÷refine; ny = _model.m_ny÷refine; 
    nelem = nx*ny; el = 0; line = ny * refine * refine - ny; buffer = 0;
    # Open and read file:
    open(_filename, "r") do io
        bin_array = read(io);        
        # Build the element material map based on the refinement level:
        for e=1:nelem
            buffer = _model.m_matKeys[bin_array[e]+1];
            for i=1:refine
                for j=1:refine
                    el = e + ((e-1)÷ny)*line + (j-1) + ((i-1)%refine)*ny*refine + ((e-1)%ny)*(refine-1);
                    _elemMatMap[el] = buffer;
                end        
            end 
        end        
        bin_array = nothing;
    end
end

# Estimate memory consuption:
"""
    estimateMemory(_model::Model)

    The amount of RAM used depends on the size of the problem and the type of solver chosen.

"""
function estimateMemory(_model::Model)
    println(".Estimate memory!");
    # elemMatMap = 16 bits * nElems
    # DOFMap = 64 bits * nNodes
    # RHS = 64 bits * nDOFs
    # PCGM Solver   / solverType == 0 / M d x q =  * 64 bits * nDOFs
    # Direct Solver / solverType == 1 / K = 18 * 64 bits * nElems (rough sparse estimative)
    mem = 0;
    if (_model.m_solverType == 0)
        mem += (16*_model.m_nElems + 64*_model.m_nNodes + 5*64*_model.m_nDOFs)/8/1_000_000;
    elseif (_model.m_solverType == 1)
        mem += (16*_model.m_nElems + 64*_model.m_nNodes + 2*64*_model.m_nDOFs + 18*64*_model.m_nElems)/8/1_000_000;
    end
    println("$(_model.m_nDOFs) DOFs");
    println("$mem MB");
end

# Compute the element conductivity matrix for each material:
"""
    elementConductivityMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})

    One of pfem4ec's differentials is the fact that it does not need to calculate all the elements of the model, 
    only the elements with different materials. This saves RAM and CPU usage.
"""
function elementConductivityMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})
    println(".Compute each element conductivity matrix!");
    # Compute the matrices for each material:
    i = 0;
    for elemProps in _model.m_matProp[_model.m_matKeys.!=0]
        i += 1;
        _K[:,:,i], _B[:,:,i] = Q4ElementConductivity(elemProps);
    end
end

# Element Q4 Conductivity - FEM:
"""
    Q4ElementConductivity(_elemProps::Float64)

    Analytical solutions for pixel-based finite element method 
    (Q4 quadrilateral element that has four nodes)
"""
function Q4ElementConductivity(_elemProps::Float64)
    # Initializations:
    K = zeros(Float64,4,4);
    B = zeros(Float64,2,4);
    # Analytical k and B for a Pixel-Based FEM: 
    a = 2/3*_elemProps; b = 1/3*_elemProps; c = 1/6*_elemProps;
    q = 1/2*_elemProps;
    K[1,1] = +a; K[1,2] = -c; K[1,3] = -b; K[1,4] = -c;
    K[2,1] = -c; K[2,2] = +a; K[2,3] = -c; K[2,4] = -b;
    K[3,1] = -b; K[3,2] = -c; K[3,3] = +a; K[3,4] = -c;
    K[4,1] = -c; K[4,2] = -b; K[4,3] = -c; K[4,4] = +a;
    B[1,1] = -q; B[1,2] = +q; B[1,3] = +q; B[1,4] = -q;
    B[2,1] = -q; B[2,2] = -q; B[2,3] = +q; B[2,4] = +q;
    return K, B
end

# Generate the Degree of Freedom Map:
"""
    generateDOFMap!(_model::Model, _DOFMap::Array{UInt64,1})

    This function create degrees of freedom for the model
    top to bottom and left to right
    The total number of nodes is (nx+1)*(ny+1)
"""
function generateDOFMap!(_model::Model, _DOFMap::Array{UInt64,1})
    println(".Generate the Map of DOFs (Degrees of Freedom)!");
    # Number the DOFs following the nodes from top to bottom and left to right:
    for n=1:_model.m_nNodes
        i = (n-1)%_model.m_nNodes;
        _DOFMap[n] = (i - (i÷(_model.m_ny+1)) - _model.m_ny*((i%(_model.m_ny+1)÷_model.m_ny)))%_model.m_nElems + 1;
    end
end

# Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1, axis 0 = X || axis 1 = Y
"""
    computeRHS!(_model::Model, _RHS::Array{Float64,1}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    There are two different boundaries conductions: Boundary or Domain

"""
function computeRHS!(_model::Model, _RHS::Array{Float64,1}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".Compute RHS!");
    # Initializations:
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; e = 0; 
    # Compute each RHS (_axis) based on boundary or domain data (_model.m_rhsType):
    if _model.m_rhsType == 1  # Boundary         
        if _axis == 0
            deltaT = _model.m_nx;
            c = _model.m_nx;
            for r=1:_model.m_ny
                e = r + (c-1)*_model.m_ny;
                N1 = e + c; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
                _RHS[_DOFMap[N1]] -= (_K[1,2,_elemMatMap[e]] + _K[1,3,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N2]] -= (_K[2,2,_elemMatMap[e]] + _K[2,3,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N3]] -= (_K[3,2,_elemMatMap[e]] + _K[3,3,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N4]] -= (_K[4,2,_elemMatMap[e]] + _K[4,3,_elemMatMap[e]])*deltaT;
            end
        elseif _axis == 1
            deltaT = _model.m_ny;
            r = 1;
            for c=1:_model.m_nx
                e = r + (c-1)*_model.m_ny;
                N1 = e + c; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
                _RHS[_DOFMap[N1]] -= (_K[1,3,_elemMatMap[e]] + _K[1,4,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N2]] -= (_K[2,3,_elemMatMap[e]] + _K[2,4,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N3]] -= (_K[3,3,_elemMatMap[e]] + _K[3,4,_elemMatMap[e]])*deltaT;
                _RHS[_DOFMap[N4]] -= (_K[4,3,_elemMatMap[e]] + _K[4,4,_elemMatMap[e]])*deltaT;
            end
        end
    elseif _model.m_rhsType == 0  # Domain      
        for e=1:_model.m_nElems
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            _RHS[_DOFMap[N1]] += _B[_axis+1,1,_elemMatMap[e]];
            _RHS[_DOFMap[N2]] += _B[_axis+1,2,_elemMatMap[e]];
            _RHS[_DOFMap[N3]] += _B[_axis+1,3,_elemMatMap[e]];
            _RHS[_DOFMap[N4]] += _B[_axis+1,4,_elemMatMap[e]];
        end    
    end
end

# Direct Solver: [K] 64 bits * _numDOFs * _numDOFs 
"""
    directMethod!(_model::Model, _x1::Array{Float64,1}, _x2::Array{Float64,1}, _RHS1::Array{Float64,1}, _RHS2::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    For the direct solver it is necessary to calculate the stiffness only once (saving the CPU usage), 
    however for large models, even in sparse format, a lot of RAM is used.

"""
function directMethod!(_model::Model, _x1::Array{Float64,1}, _x2::Array{Float64,1}, _RHS1::Array{Float64,1}, _RHS2::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".Direct Solver!");
    # Initializations:
    K = spzeros(_model.m_nDOFs,_model.m_nDOFs);
    pElemDOFNum = zeros(UInt64,4);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    # Assembly system matrix:
    for e=1:_model.m_nElems
        N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
        pElemDOFNum[1] = _DOFMap[N1]; pElemDOFNum[2] = _DOFMap[N2]; pElemDOFNum[3] = _DOFMap[N3]; pElemDOFNum[4] = _DOFMap[N4];
        for i=1:4
            for j=1:4
                K[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j,_elemMatMap[e]];
            end
        end
    end
    # Solve for two rhs:
    _x1 .= K\_RHS1;
    _x2 .= K\_RHS2;
end

# Jacobi Preconditioner: assembly || M
"""
    jacobiPrecond!(_model::Model, _M::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    The Conjugated Gradients method in its direct form, however, presents problems of numerical instability.
    To alleviate the problem, the preconditioner M is used.
    Jacobi's preconditioner is the diagonal matrix of stiffness. 

"""
function jacobiPrecond!(_model::Model, _M::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".Jacobi Preconditioner!");
    # Initializations:
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    # Compute the preconditioner: 
    for e=1:_model.m_nElems
        N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
        _M[_DOFMap[N1]] += _K[1,1,_elemMatMap[e]];
        _M[_DOFMap[N2]] += _K[2,2,_elemMatMap[e]];
        _M[_DOFMap[N3]] += _K[3,3,_elemMatMap[e]];
        _M[_DOFMap[N4]] += _K[4,4,_elemMatMap[e]];
    end
    _M .= _M.\1;
end 

# Preconditioned Conjugate Gradient Method:
"""
    pcg!(_model::Model, _x::Array{Float64,1}, _r::Array{Float64,1}, _M::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    The conjugate gradient method is an algorithm for the numerical solution of
    particular systems of linear equations, those whose matrix is defined symmetric and positive.

    This Preconditioned Conjugate Gradient method is based in [Shewchuk algorithm](https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf)

"""
function pcg!(_model::Model, _x::Array{Float64,1}, _r::Array{Float64,1}, _M::Array{Float64,1}, _K::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".PCG Solver!");
    # Initializations:    
    d = zeros(Float64,_model.m_nDOFs);
    q = zeros(Float64,_model.m_nDOFs);
    pElemDOFNum = zeros(UInt64,4);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    q_temp = 0;    
    # PCG Initialization:
    d .= _r;
    d .*= _M;
    delta_new = (_r'*d)[1,1];
    delta_0 = delta_new;
    i_max = _model.m_pcgIter;
    ii = 0;
    # PCG Iterations:
    while (ii<i_max) && (abs(delta_new)>_model.m_pcgTol*_model.m_pcgTol*abs(delta_0)) #(maximum(abs.(_r))>_pcgTol)
        @fastmath @inbounds @simd for e=1:_model.m_nElems
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = _DOFMap[N1]; pElemDOFNum[2] = _DOFMap[N2]; pElemDOFNum[3] = _DOFMap[N3]; pElemDOFNum[4] = _DOFMap[N4];
            for i=1:4
                q_temp = 0;
                for j=1:4
                    q_temp += _K[i,j,_elemMatMap[e]] * d[pElemDOFNum[j]];
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        alfa = delta_new/(d'*q)[1,1];
        d .*= alfa;
        _x .+= d;
        q .*= alfa;
        _r .-= q;
        q .= _r;
        q .*= _M;
        delta_old = delta_new;
        delta_new = (_r'*q)[1,1];
        beta = delta_new/delta_old;
        d .*= beta/alfa;
        d .+= q;
        q .*= 0;
        ii += 1;
    end
    println("$ii steps");
    println(sqrt(abs(delta_new)/abs(delta_0)));
end

# Compute Flux-FEM Effective property:
"""
    femEffective(_model::Model, _T::Array{Float64,1}, _axis::Int, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    After solving the system of equations, this function calculates the effective property of the material.

"""
function femEffective(_model::Model, _T::Array{Float64,1}, _axis::Int, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".Compute Effective Property!");
    # Initializations:
    QX = 0; QY = 0;
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    pElemDOFNum = zeros(UInt64,4);
    C = zeros(Float64,2,2);
    # Compute the effective properties for each test: 
    if _model.m_rhsType == 1  # Boundary
        if _axis == 0
            deltaT = _model.m_nx;   
            for eb = _model.m_nElems-(_model.m_ny-1):_model.m_nElems  
                QX += (_B[1,2,_elemMatMap[eb]]*deltaT); QX += (_B[1,3,_elemMatMap[eb]]*deltaT);
                QY += (_B[2,2,_elemMatMap[eb]]*deltaT); QY += (_B[2,3,_elemMatMap[eb]]*deltaT);             
            end    
        elseif _axis == 1
            deltaT = _model.m_ny;
            for eb = 1:(_model.m_ny):_model.m_nElems
                QX += (_B[1,3,_elemMatMap[eb]]*deltaT); QX += (_B[1,4,_elemMatMap[eb]]*deltaT);
                QY += (_B[2,3,_elemMatMap[eb]]*deltaT); QY += (_B[2,4,_elemMatMap[eb]]*deltaT);             
            end 
        end
        for e=1:_model.m_nElems
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i=1:4
                QX += (_B[1,i,_elemMatMap[e]]*_T[_DOFMap[pElemDOFNum[i]]]);
                QY += (_B[2,i,_elemMatMap[e]]*_T[_DOFMap[pElemDOFNum[i]]]);
            end
        end        
    elseif _model.m_rhsType == 0  # Domain
        t = zeros(Float64,4);
        if (_axis == 0);     t[1] = 0; t[2] = 1; t[3] = 1; t[4] = 0;
        elseif (_axis == 1); t[1] = 0; t[2] = 0; t[3] = 1; t[4] = 1; end   
        for e = 1:_model.m_nElems            
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i=1:4
                QX += (_B[1,i,_elemMatMap[e]]*(t[i] -_T[_DOFMap[pElemDOFNum[i]]]));
                QY += (_B[2,i,_elemMatMap[e]]*(t[i] -_T[_DOFMap[pElemDOFNum[i]]]));
            end
        end 
    end
    C[1,_axis+1] = QX/_model.m_nElems; C[2,_axis+1] = QY/_model.m_nElems;
    return C
end

# Compute Flux-FEM Effective property:
"""
    femEffective(_model::Model, _Tx::Array{Float64,1}, _Ty::Array{Float64,1}, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})

    After solving the system of equations, this function calculates the effective property of the material.

"""
function femEffective(_model::Model, _Tx::Array{Float64,1}, _Ty::Array{Float64,1}, _B::Array{Float64,3}, _DOFMap::Array{UInt64,1}, _elemMatMap::Array{UInt16,1})
    println(".Compute Effective Property!");
    # Initializations:
    QXx = 0; QXy  = 0; QYx = 0; QYy = 0;
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    pElemDOFNum = zeros(UInt64,4);
    # Compute the effective properties for all tests:
    if _model.m_rhsType == 1  # Boundary
        deltaTx = _model.m_nx;   
        for eb=_model.m_nElems-(_model.m_ny-1):_model.m_nElems  
            QXx += (_B[1,2,_elemMatMap[eb]]*deltaTx); QXx += (_B[1,3,_elemMatMap[eb]]*deltaTx);
            QXy += (_B[2,2,_elemMatMap[eb]]*deltaTx); QXy += (_B[2,3,_elemMatMap[eb]]*deltaTx);          
        end
        deltaTy = _model.m_ny;
        for eb=1:(_model.m_ny):_model.m_nElems
            QYx += (_B[1,3,_elemMatMap[eb]]*deltaTy); QYx += (_B[1,4,_elemMatMap[eb]]*deltaTy);
            QYy += (_B[2,3,_elemMatMap[eb]]*deltaTy); QYy += (_B[2,4,_elemMatMap[eb]]*deltaTy);             
        end 
        for e=1:_model.m_nElems
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i=1:4
                QXx += (_B[1,i,_elemMatMap[e]]*_Tx[_DOFMap[pElemDOFNum[i]]]);
                QXy += (_B[2,i,_elemMatMap[e]]*_Tx[_DOFMap[pElemDOFNum[i]]]);
                QYx += (_B[1,i,_elemMatMap[e]]*_Ty[_DOFMap[pElemDOFNum[i]]]);
                QYy += (_B[2,i,_elemMatMap[e]]*_Ty[_DOFMap[pElemDOFNum[i]]]);
            end
        end
    elseif _model.m_rhsType == 0  # Domain
        t1 = zeros(Float64,4); t2 = zeros(Float64,4);
        t1[1] = 0; t1[2] = 1; t1[3] = 1; t1[4] = 0;
        t2[1] = 0; t2[2] = 0; t2[3] = 1; t2[4] = 1;        
        for e=1:_model.m_nElems
            N1 = e + ((e-1)÷_model.m_ny) + 1; N3 = N1 + _model.m_ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i=1:4
                QXx += (_B[1,i,_elemMatMap[e]]*(t1[i] -_Tx[_DOFMap[pElemDOFNum[i]]]));
                QXy += (_B[2,i,_elemMatMap[e]]*(t1[i] -_Tx[_DOFMap[pElemDOFNum[i]]]));
                QYx += (_B[1,i,_elemMatMap[e]]*(t2[i] -_Ty[_DOFMap[pElemDOFNum[i]]]));
                QYy += (_B[2,i,_elemMatMap[e]]*(t2[i] -_Ty[_DOFMap[pElemDOFNum[i]]]));
            end
        end
    end
    C = [ QXx/_model.m_nElems QYx/_model.m_nElems; QXy/_model.m_nElems QYy/_model.m_nElems];
    return C
end

# -----------------
# Main function
"""
    pfem4ec(_arg)

    The main function calls all the other functions:
    read files, create model, generate boundary condictions,
    solve problem, compute effective properties
    
"""
function pfem4ec(_arg)
    # Read JSON File and build the Model data structure: 
    model = readJSON!(_arg*".json");
    # Read RAW File and setup the Element Material Map:
    elemMatMap = zeros(UInt16,model.m_nElems);
    readRAW!(model,elemMatMap,_arg*".raw");
    # Estimate Memory Consumption:
    estimateMemory(model);
    # Compute the conductivity matrix for each Material:
    mK = zeros(Float64,4,4,model.m_nMat);
    mB = zeros(Float64,2,4,model.m_nMat);
    elementConductivityMatrices!(model,mK,mB);
    # Generate the Degree of Freedom Map:
    vDOFMap = zeros(UInt64,model.m_nNodes); 
    generateDOFMap!(model,vDOFMap);
    # SOLVE:
    if (model.m_solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:
        C = zeros(Float64,2,2);
        vRHS = zeros(Float64,model.m_nDOFs);   
        vX = zeros(Float64,model.m_nDOFs);
        vM = zeros(Float64,model.m_nDOFs);
        # Compute the Jacobi preconditioner:
        jacobiPrecond!(model,vM,mK,vDOFMap,elemMatMap);
        for axis=0:1
            # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
            computeRHS!(model,vRHS,axis,mK,mB,vDOFMap,elemMatMap);            
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):    
            GC.gc();
            pcg!(model,vX,vRHS,vM,mK,vDOFMap,elemMatMap);
            GC.gc();        
            # Compute Effective Property:
            C .+= femEffective(model,vX,axis,mB,vDOFMap,elemMatMap);
            vRHS .*= 0;
            vX .*= 0;
        end
        M = nothing;
        println(C)
    elseif (model.m_solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
        vRHS1 = zeros(Float64,model.m_nDOFs); vRHS2 = zeros(Float64,model.m_nDOFs);
        computeRHS!(model,vRHS1,0,mK,mB,vDOFMap,elemMatMap);
        computeRHS!(model,vRHS2,1,mK,mB,vDOFMap,elemMatMap);
        # Solver
        vX1 = zeros(Float64,model.m_nDOFs); vX2 = zeros(Float64,model.m_nDOFs);
        directMethod!(model,vX1,vX2,vRHS1,vRHS2,mK,vDOFMap,elemMatMap);
        vRHS1 = nothing; vRHS2 = nothing;
        # Compute Effective Property:
        C = femEffective(model,vX1,vX2,mB,vDOFMap,elemMatMap);
        vX1 = nothing; vX2 = nothing;
        println(C)
    end
    println("--------------------------------------");
    vDOFMap = nothing; elemMatMap = nothing;
end

# Starts application
if length(ARGS) > 0
    @time pfem4ec(ARGS[1])
end
