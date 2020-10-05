# Compute Effective Properties

using LinearAlgebra
using JSON
using SparseArrays


# Functions
# -----------------------------------------------------

# -----------------------------------------------------
# Read JSON file
function readJSON(filename::String)

    println("----")
    println("Read JSON")
    
    # Initialize
    m_materials_I = zeros(UInt16,256);
    m_materials_F = zeros(Float64,256);
    
    data = Dict();
    open(filename,"r") do f
        data = JSON.parse(f);
    end
    
    m_tol = data["TOL"];
    m_nx = data["W"];
    m_ny = data["H"];
    m_solver = data["SOLVER"];
    c_RHS = data["RHS"];
    m_numMat = data["NMAT"];
    mat = data["MATERIALS"];
    m_nrefine = data["REFINE"];
    for i=1:m_numMat
        aux = mat[i];
        m_materials_I[convert(UInt8,aux[1])+1] = i;
        m_materials_F[convert(UInt8,aux[1])+1] = convert(Float64, aux[2]);
    end
    
    return m_tol, m_materials_I, m_materials_F, m_numMat, m_nx, m_ny, m_nrefine, m_solver, c_RHS
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Read RAW file
    function readRAW(filename::String, m_nx::Int, m_ny::Int, m_nrefine::Int)
    
    println("----")
    println("Read RAW")
    
    io = open(filename, "r");
    bin_array = read(io);
    
    m_matID = reshape(reinterpret(UInt8, bin_array),:);
    m_matID = reshape(m_matID,m_nx,m_ny);
    m_matID = convert(Array{UInt64},kron(m_matID,ones(m_nrefine,m_nrefine))).+1;
    m_matID = m_matID[:];
    m_nx = m_nx*m_nrefine; 
    m_ny = m_nx;
    
    return m_matID, m_nx, m_ny
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Conductivity Matrix Materials
    function matsCondMatrix(m_gdlNo::Int,m_numMat::Int,m_materials_I::Array{UInt16,1},m_materials_F::Array{Float64,1})
    
    println("----")
    println("Initialize Variables")
    
    m_k = zeros(m_gdlNo*4,m_gdlNo*4,m_numMat);
    m_B = zeros(m_gdlNo*2,m_gdlNo*4,m_numMat);
    
    i = 0;
    for prop in m_materials_F[m_materials_I.!=0]
        i += 1;
        elems_prop = prop;
        m_k[:,:,i] , m_B[:,:,i] = Q4ElementConductivity(elems_prop);
    end
    
    return m_k, m_B
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Element Q4 Conductivity - FEM
    function Q4ElementConductivity(elems_prop::Float64)
    
    # Element
    x = [0.;1.;1.;0.];
    y = [0.;0.;1.;1.];
    
    # Constitutive
    c = elems_prop;
    Ident = Matrix{Float64}(I, 2, 2);
    C = c*Ident;
    
    # Gauss Points
    PG = [-1.0/sqrt(3) 1.0/sqrt(3)];
    w = [1.0 1.0];
    
    k = zeros(4,4);
    BC = zeros(2,4);
    
    for i=1:2
        r = PG[1,i];
        wx = w[1,i];
        for j = 1:2
            s = PG[1,j];
            wy = w[1,j];
            B,J = Q4BMatrix(r, s, x, y);
            dJ = det(J);
            k  += B'*C*B*dJ*wx*wy;
            BC += C*B*dJ*wx*wy;
        end
    end
    
    return k, BC
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Q4BMatrix - FEM
    function Q4BMatrix(r::Float64, s::Float64, x::Array{Float64,1}, y::Array{Float64,1})
    
    X = [x'; y'];
    
    dN1dr = -(1-r)*.25;
    dN2dr =  (1-r)*.25;
    dN3dr =  (1+r)*.25;
    dN4dr = -(1+r)*.25;
    
    dN1ds = -(1-s)*.25;
    dN2ds = -(1+s)*.25;
    dN3ds =  (1+s)*.25;
    dN4ds =  (1-s)*.25;
    
    dN = [dN1dr dN2dr dN3dr dN4dr;
          dN1ds dN2ds dN3ds dN4ds];
    
    J = dN*X';
    
    dNdx = J\dN;
    
    B = [dNdx[1,1] dNdx[1,2] dNdx[1,3] dNdx[1,4];
         dNdx[2,1] dNdx[2,2] dNdx[2,3] dNdx[2,4]];
    
    return B, J
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Degree of Freedom Map
    function get_dofmap(m_nx::Int, m_ny::Int, m_numElem::Int)
    
    println("----")
    println("Degree of Freedom Map")
    
    dofMap = zeros(UInt64,m_ny+1,m_nx+1);
    
    dofMap[1:m_ny,1:m_nx] = reshape(1:m_numElem,m_ny,m_nx);
    dofMap[end,:] = dofMap[1,:];
    dofMap[:,end] = dofMap[:,1];
    dofMap = dofMap[:];
    
    return dofMap
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Compute RHS
    # c_RHS Boundary = 0 || Domain = 1
    # axis x = 0 || y = 1
    function computeRHS(dofMap::Array{UInt64,1}, c_RHS::Int, axis::Int, m_k::Array{Float64,3}, m_nGDL::Int, m_nx::Int, m_ny::Int, m_numElem::Int, m_matID::Array{UInt64,1} ,m_materials_I::Array{UInt16,1})
    
    println("----")
    println("Compute RHS")
    m_RHS = zeros(m_nGDL);
    
    if c_RHS == 0     # Boundary
    
        if axis == 0
    
            deltaT = m_nx;
            c = m_nx;
            for r=1:m_ny
                elem = r + (c-1)*m_ny;
                k = m_k[:,:,m_materials_I[m_matID[elem]]];
                N1 = r + (c-1)*m_ny + c;
                N2 = r + c*m_ny + c + 1;
                N3 = r + c*m_ny + c;
                N4 = r + (c-1)*m_ny + c - 1;
    
                m_RHS[dofMap[N1]] += - (k[1,2] + k[1,3])*deltaT;
                m_RHS[dofMap[N2]] += - (k[2,2] + k[2,3])*deltaT;
                m_RHS[dofMap[N3]] += - (k[3,2] + k[3,3])*deltaT;
                m_RHS[dofMap[N4]] += - (k[4,2] + k[4,3])*deltaT;
            end
    
        elseif axis == 1
    
            deltaT = m_ny;
            r = 1;
            for c=1:m_nx
                elem = r + (c-1)*m_ny;
                k = m_k[:,:,m_materials_I[m_matID[elem]]];
                N1 = r + (c-1)*m_ny + c;
                N2 = r + c*m_ny + c + 1;
                N3 = r + c*m_ny + c;
                N4 = r + (c-1)*m_ny + c - 1;
    
                m_RHS[dofMap[N1]] += - (k[1,3] + k[1,4])*deltaT;
                m_RHS[dofMap[N2]] += - (k[2,3] + k[2,4])*deltaT;
                m_RHS[dofMap[N3]] += - (k[3,3] + k[3,4])*deltaT;
                m_RHS[dofMap[N4]] += - (k[4,3] + k[4,4])*deltaT;
            end
    
        end
    
    elseif c_RHS == 1  # Domain
    
        for e=1:m_numElem
            c = floor(UInt64, ((e-1)/m_ny) + 1);
            b = m_B[:,:,m_materials_I[m_matID[e]]];
            N1 = e + c;
            N2 = e + c + 1 + m_ny;
            N3 = e + c + m_ny;
            N4 = e + c - 1;
    
            if axis == 0
                m_RHS[m_dofMap[N1]] += b[1,1];
                m_RHS[m_dofMap[N2]] += b[1,2];
                m_RHS[m_dofMap[N3]] += b[1,3];
                m_RHS[m_dofMap[N4]] += b[1,4];
            elseif axis == 1
                m_RHS[m_dofMap[N1]] += b[2,1];
                m_RHS[m_dofMap[N2]] += b[2,2];
                m_RHS[m_dofMap[N3]] += b[2,3];
                m_RHS[m_dofMap[N4]] += b[2,4];
            end
        end
    end
    
    return m_RHS
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Pre Conjugate Gradient Solver
    function pcg(dofMap::Array{UInt64,1}, r::Array{Float64,1}, m_k::Array{Float64,3}, m_tol::Float64, m_nGDL::Int, m_nx::Int, m_ny::Int, m_numElem::Int, m_matID::Array{UInt64,1}, m_materials_I::Array{UInt16,1})
    
    println("----")
    println("PCG Solver")
    
    # Preconditioner Jacobi assembly || M
    # -----------------------------------------------------
    
    M = zeros(m_nGDL);
    x = zeros(m_nGDL);
    edof1 = 4;
    edof2 = 4;
    
    for e=1:m_numElem
        c = floor(UInt64, ((e-1)/m_ny) + 1);
        k = m_k[:,:,m_materials_I[m_matID[e]]];
        N1 = e + c;
        N2 = e + c + 1 + m_ny;
        N3 = e + c + m_ny;
        N4 = e + c - 1;
        M[dofMap[N1]] += k[1,1];
        M[dofMap[N2]] += k[2,2];
        M[dofMap[N3]] += k[3,3];
        M[dofMap[N4]] += k[4,4];
    end
    
    # -----------------------------------------------------
    # Compute r
    # -----------------------------------------------------
    
    # loop in elements to compute r
    for e=1:m_numElem
        c = floor(UInt64, ((e-1)/m_ny) + 1);
        k = m_k[:,:,m_materials_I[m_matID[e]]];
    
        N1 = e + c;
        N2 = e + c + 1 + m_ny;
        N3 = e + c + m_ny;
        N4 = e + c - 1;
        pElemDOFNum = [dofMap[N1] dofMap[N2] dofMap[N3] dofMap[N4]];
    
        for i=1:edof1
            r_temp = 0;
            for j=1:edof2
                r_temp += k[i,j] * x[pElemDOFNum[j]];
            end
            r[pElemDOFNum[i]] -= r_temp;
        end
    end
    
    # -----------------------------------------------------
    # Preconditioned Conjugate Gradient
    # -----------------------------------------------------
    
    # initialize
    d = M.\r;
    delta_new = transpose(r)*d;
    # delta_0 = delta_new
    ii = 0;
    i_max = 2000;
    
    # while (ii<i_max) && (abs(delta_new)>eps*eps*abs(delta_0))
    while (ii<=i_max) && (maximum(abs.(r))>m_tol)
    
        q = zeros(m_nGDL);
    
        for e=1:m_numElem
            c = floor(UInt64, ((e-1)/m_ny) + 1);
            k = m_k[:,:,m_materials_I[m_matID[e]]];
    
            N1 = e + c;
            N2 = e + c + 1 + m_ny;
            N3 = e + c + m_ny;
            N4 = e + c - 1;
            pElemDOFNum = [dofMap[N1] dofMap[N2] dofMap[N3] dofMap[N4]];
    
            for i=1:edof1
                q_temp = 0;
                for j=1:edof2
                    q_temp += k[i,j] * d[pElemDOFNum[j]];
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
    
        q_temp = nothing;
    
        alfa = delta_new/(transpose(d)*q);
        x += d*alfa;
        #    if rem(ii,50) == 0
        #       r = b - A*x;
        #    else
        #       r = r - alfa*q;
        #    end
        r -= q*alfa;
        q = nothing;
        s = M.\r;
        delta_old = delta_new;
        delta_new = transpose(r)*s;
        beta = delta_new/delta_old;
        d = s + d*beta;
        ii += 1;
        s = nothing;
    
    end
    
    return x
    
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Recover Temperature Values
    # axis x = 0 || y = 1
    function rcvTemp( T1::Array{Float64,1}, c_RHS::Int, axis::Int, m_nx::Int, m_ny::Int, m_numNos::Int)
    
    println("----")
    println("Recover Temperature Values")
    
    T0 = zeros(m_numNos);
    
    if c_RHS == 0    # Boundary
    
        if axis == 0
            deltaT = m_nx;
            T0[(m_numNos-m_ny):m_numNos] .= deltaT;  
        elseif axis == 1
            deltaT = m_ny;
            T0[1:(m_ny+1):m_numNos] .= deltaT;     
        end
        T = T1 + T0;
    
    elseif c_RHS == 1 # Domain
    
        if axis == 0
            flux = 0:m_nx;
            T0 = vcat(fill.(flux, (m_ny+1))...);
        elseif axis == 1
            flux = m_ny:-1:0;
            T0 = vcat(vcat(fill.(flux, (m_ny+1))'...)...);
        end
        T = T0 - T1;
    end
    
    return T
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    # Compute Flux - FEM
    # Effective property
    function femEffective(Tx::Array{Float64,1},Ty::Array{Float64,1}, m_B::Array{Float64,3}, m_nx::Int, m_ny::Int, m_numElem::Int, m_matID::Array{UInt64,1}, m_materials_I::Array{UInt16,1})
    
    println("----")
    println("Compute Effective Property")
    
    QXx = 0; QXy  = 0; QYx = 0; QYy = 0 ; V = 0;
    for e=1:m_numElem
        c = floor(UInt64, ((e-1)/m_ny) + 1);
        N1 = e + c;
        N2 = e + c + 1 + m_ny;
        N3 = e + c + m_ny;
        N4 = e + c - 1;
    
        pElemDOFNum = [N1 N2 N3 N4];
        b = m_B[:,:,m_materials_I[m_matID[e],1],1];
    
        q1 = b*Tx[pElemDOFNum]';
        q2 = b*Ty[pElemDOFNum]';
        QXx = QXx + q1[1];
        QXy = QXy + q1[2];
        QYx = QYx + q2[1];
        QYy = QYy + q2[2];
        V = V + 1;
    end
    
    C = [ QXx/V QXy/V; QYx/V QYy/V];
    
    return C
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    function memEstimate(m_numNos::Int, m_numElem::Int, m_nGDL::Int, m_solver::Int)
        
    println("----")
    println("Memory Estimate")
    # m_matID = 64 bits * m_numElem
    # m_dofMap = 64 bits * m_numNos
    # m_RHS, m_RHS_2 = 2 * 64 bits * m_nGDL
    # PCG - m_solver == 0
    # M r d x q = 5 * 64 bits * m_nGDL
    # Direct Solver - m_solver == 1
    # K = 64 bits * 16 * m_numElem (rough sparse estimative)
    
    mem = 125.0;
    if (m_solver == 0)
        mem += (64*m_numElem + 64*m_numNos + 2*64*m_nGDL + 6*64*m_nGDL)/8/1000/1000;
    elseif (m_solver == 1)
        mem += (64*m_numElem + 64*m_numNos + 2*64*m_nGDL + 18*64*m_numElem)/8/1000/1000;
    end
    @printf "%.2f MB\n" mem
    
    end
    # -----------------------------------------------------
    
    # -----------------------------------------------------
    
        
    # -----------------
    # Read File | JSON
    # -----------------
    filename = ARGS[1] * ".json";
    m_tol , m_materials_I, m_materials_F, m_numMat, m_nx, m_ny, m_nrefine, m_solver, c_RHS = readJSON(filename);
    #println(varinfo())
    
    # ----------------
    # Read File | RAW
    # ----------------
    filename = ARGS[1] * ".raw";
    m_matID , m_nx, m_ny = readRAW(filename,m_nx,m_ny,m_nrefine);
    #println(varinfo())
    
    # ------------------
    # Initial Variables
    # ------------------
     m_gdlNo = 1;
     m_numNos = (m_nx+1)*(m_ny+1);
     m_numElem = (m_nx)*(m_ny);
     m_nGDL = m_gdlNo*m_nx*m_ny;
    
     # ----------------
     # Memory Estimate
     # ----------------
     memEstimate(m_numNos, m_numElem, m_nGDL, m_solver);
    
     # ------------------------------
     # Conductivity Matrix Materials
     # ------------------------------
     m_k, m_B = matsCondMatrix(m_gdlNo, m_numMat, m_materials_I, m_materials_F);
    #println(varinfo())
    
     # ----------------------
     # Degree of Freedom Map
     # ----------------------
     m_dofMap = get_dofmap(m_nx, m_ny, m_numElem);
    #println(varinfo())
    
     # ---------------------------------
     # Compute RHS - Boundary or Domain
     # ---------------------------------
     # c_RHS Boundary = 0 || Domain = 1
     # axis 0 = X || axis 1 = Y
     m_RHS   = computeRHS(m_dofMap, c_RHS, 0, m_k, m_nGDL, m_nx, m_ny, m_numElem, m_matID, m_materials_I);
     m_RHS_2 = computeRHS(m_dofMap, c_RHS, 1, m_k, m_nGDL, m_nx, m_ny, m_numElem, m_matID, m_materials_I);
    
     # --------
     # SOLVERS
     # --------
    println("----")
    println("Direct Solver")
        
    K = spzeros(m_nGDL,m_nGDL);
    #println(varinfo())
    edof1 = 4; edof2 = 4;
    for e=1:m_numElem
        c = floor(UInt64, ((e-1)/m_ny) + 1);
        k = m_k[:,:,m_materials_I[m_matID[e],1],1];
    
        N1 = e + c;
        N2 = e + c + 1 + m_ny;
        N3 = e + c + m_ny;
        N4 = e + c - 1;
        pElemDOFNum = [m_dofMap[N1] m_dofMap[N2] m_dofMap[N3] m_dofMap[N4]];
    
        for i=1:edof1
            for j=1:edof2
                K[pElemDOFNum[i],pElemDOFNum[j]] += k[i,j];
            end
        end
    end
    
    tx = K\m_RHS;
    ty = K\m_RHS_2;
     m_RHS = nothing; m_RHS_2 = nothing;
    println(varinfo())
    
     # ---------------------------
     # Recover Temperature Values
     # ---------------------------
     Tx = rcvTemp(tx[m_dofMap], c_RHS, 0, m_nx, m_ny, m_numNos);
     Ty = rcvTemp(ty[m_dofMap], c_RHS, 1, m_nx, m_ny, m_numNos);
     tx = nothing; ty = nothing; m_dofMap = nothing;
    
     # ---------------------------
     # Compute Effective Property
     # ---------------------------
     C = femEffective(Tx,Ty,m_B,m_nx,m_ny,m_numElem,m_matID,m_materials_I);
     println(C)
    
    
    
    
    