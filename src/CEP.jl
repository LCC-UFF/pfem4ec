# Compute Effective Properties
using LinearAlgebra

# Functions
# -----------------------------------------------------

# -----------------------------------------------------
# Element Q4 Conductivity - FEM
function Q4ElementConductivity(elems_prop)

    # Element
    x = [0;1;1;0]
    y = [0;0;1;1]
    
    # Constitutive
    c = elems_prop
    Ident = Matrix{Float64}(I, 2, 2)
    C = c*Ident
    
    # Gauss Points
    PG = [-1/sqrt(3) 1/sqrt(3)]
    w = [1 1]
    
    k = zeros(4,4)
    BC = zeros(2,4)
    
    for i=1:2
        r = PG[1,i]
        wx = w[1,i]
        for j = 1:2
            s = PG[1,j]
            wy = w[1,j]
            B,J = Q4BMatrix(r, s, x, y)
            dJ = det(J)
            k  += B'*C*B*dJ*wx*wy
            BC += C*B*dJ*wx*wy;
        end
    end

return k, BC
end
# -----------------------------------------------------
    
# -----------------------------------------------------
# Q4BMatrix - FEM
function Q4BMatrix(r, s, x, y)
    
    X = [x'; y']
    
    dN1dr = -(1-r)*.25
    dN2dr =  (1-r)*.25
    dN3dr =  (1+r)*.25
    dN4dr = -(1+r)*.25
    
    dN1ds = -(1-s)*.25
    dN2ds = -(1+s)*.25
    dN3ds =  (1+s)*.25
    dN4ds =  (1-s)*.25
    
    dN = [dN1dr dN2dr dN3dr dN4dr;
          dN1ds dN2ds dN3ds dN4ds]
    
    J = dN*X'
    
    dNdx = J\dN
    
    B = [dNdx[1,1] dNdx[1,2] dNdx[1,3] dNdx[1,4];
         dNdx[2,1] dNdx[2,2] dNdx[2,3] dNdx[2,4]]
    
return B, J
end
# -----------------------------------------------------

# -----------------------------------------------------
# Degree of Freedom Map
function dof_modes()

dofMap = zeros(Int64,m_ny+1,m_nx+1)

dofMap[1:m_ny,1:m_nx] = reshape(1:m_numElem,m_ny,m_nx)
dofMap[end,:] = dofMap[1,:]  
dofMap[:,end] = dofMap[:,1]
dofMap = dofMap[:]

return dofMap
end
# -----------------------------------------------------

# -----------------------------------------------------
# Compute RHS
# c_RHS Boundary = 0 || Domain = 1
# axis x = 0 || y = 1
function computeRHS(dofMap,c_RHS,axis)
    
    m_RHS = zeros(m_nGDL,1)
    
    if c_RHS == 0     # Boundary 
        
        if axis == 0      
            
            deltaT = m_nx
            c = m_nx
            for r=1:m_ny
                elem = r + (c-1)*m_ny
                k = m_k[:,:,m_materials[m_matID[elem],1],1]
                N1 = r + (c-1)*m_ny + c
                N2 = r + c*m_ny + c + 1
                N3 = r + c*m_ny + c
                N4 = r + (c-1)*m_ny + c - 1
                
                m_RHS[dofMap[N1],1] += - (k[1,2] + k[1,3])*deltaT
                m_RHS[dofMap[N2],1] += - (k[2,2] + k[2,3])*deltaT
                m_RHS[dofMap[N3],1] += - (k[3,2] + k[3,3])*deltaT
                m_RHS[dofMap[N4],1] += - (k[4,2] + k[4,3])*deltaT   
            end
                    
        elseif axis == 1
            
            deltaT = m_ny
            r = 1
            for c=1:m_nx
                elem = r + (c-1)*m_ny
                k = m_k[:,:,m_materials[m_matID[elem],1],1]
                N1 = r + (c-1)*m_ny + c
                N2 = r + c*m_ny + c + 1
                N3 = r + c*m_ny + c
                N4 = r + (c-1)*m_ny + c - 1
            
                m_RHS[dofMap[N1],1] += - (k[1,3] + k[1,4])*deltaT
                m_RHS[dofMap[N2],1] += - (k[2,3] + k[2,4])*deltaT
                m_RHS[dofMap[N3],1] += - (k[3,3] + k[3,4])*deltaT
                m_RHS[dofMap[N4],1] += - (k[4,3] + k[4,4])*deltaT  
            end    

        end
        
    elseif c_RHS == 1  # Domain
       
        for e=1:m_numElem
            c = floor(Int, ((e-1)/m_ny) + 1)
            b = m_B[:,:,m_materials[m_matID[e],1],1]
            N1 = e + c
            N2 = e + c + 1 + m_ny
            N3 = e + c + m_ny
            N4 = e + c - 1
            
            if axis == 0  
                m_RHS[m_dofMap[N1],1] += b[1,1]
                m_RHS[m_dofMap[N2],1] += b[1,2]
                m_RHS[m_dofMap[N3],1] += b[1,3]
                m_RHS[m_dofMap[N4],1] += b[1,4]
            elseif axis == 1
                m_RHS[m_dofMap[N1],1] += b[2,1]
                m_RHS[m_dofMap[N2],1] += b[2,2]
                m_RHS[m_dofMap[N3],1] += b[2,3]
                m_RHS[m_dofMap[N4],1] += b[2,4]
            end
        end
                
    end

return m_RHS
end
# -----------------------------------------------------

# -----------------------------------------------------    
function pcg(dofMap, r)
   
    
# Preconditioner Jacobi assembly || M
# -----------------------------------------------------
    
    M = zeros(m_nGDL,1)

        for e=1:m_numElem
            c = floor(Int, ((e-1)/m_ny) + 1)
            k = m_k[:,:,m_materials[m_matID[e],1],1]
            N1 = e + c
            N2 = e + c + 1 + m_ny
            N3 = e + c + m_ny
            N4 = e + c - 1
            M[dofMap[N1],1] += k[1,1]
            M[dofMap[N2],1] += k[2,2]
            M[dofMap[N3],1] += k[3,3]
            M[dofMap[N4],1] += k[4,4]
        end

# -----------------------------------------------------
# Compute r
# -----------------------------------------------------
    x = zeros(m_nGDL,1)
    edof1 = 4
    edof2 = 4
    
        # loop in elements to compute r
        for e=1:m_numElem
            c = floor(Int, ((e-1)/m_ny) + 1)
            k = m_k[:,:,m_materials[m_matID[e],1],1]
            
            N1 = e + c
            N2 = e + c + 1 + m_ny
            N3 = e + c + m_ny
            N4 = e + c - 1
            pElemDOFNum = [dofMap[N1] dofMap[N2] dofMap[N3] dofMap[N4]]
            
            for i=1:edof1
                r_temp = 0;
                for j=1:edof2
                    r_temp += k[i,j] * x[pElemDOFNum[j]]
                end
                r[pElemDOFNum[i]] -= r_temp
            end
        end
        
# -----------------------------------------------------
# Preconditioned Conjugate Gradient
# -----------------------------------------------------

    # initialize
    d = M.\r
    delta_new = transpose(r)*d
    # delta_0 = delta_new
    ii = 0
    i_max = 2000
    
    # while (ii<i_max) && (abs(delta_new)>eps*eps*abs(delta_0))
    while (ii<=i_max) && (maximum(abs.(r))>m_tol)
        
        q = zeros(m_nGDL,1)
                          
        for e=1:m_numElem
            c = floor(Int, ((e-1)/m_ny) + 1)
            k = m_k[:,:,m_materials[m_matID[e],1],1]
                
            N1 = e + c
            N2 = e + c + 1 + m_ny
            N3 = e + c + m_ny
            N4 = e + c - 1
            pElemDOFNum = [dofMap[N1] dofMap[N2] dofMap[N3] dofMap[N4]]
                
            for i=1:edof1
                q_temp = 0;
                for j=1:edof2
                    q_temp += k[i,j] * d[pElemDOFNum[j]];                        
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        
        alfa = delta_new/(transpose(d)*q);
        x = x + d*alfa;
        #    if rem(ii,50) == 0
        #       r = b - A*x;
        #    else
        #       r = r - alfa*q;
        #    end
        r = r - q*alfa;
        s = M.\r;
        delta_old = delta_new;
        delta_new = transpose(r)*s;
        beta = delta_new/delta_old;
        d = s + d*beta;    
        ii += 1;
        
    end

return x

end
# -----------------------------------------------------

# -----------------------------------------------------
# Recover Temperature Values
# axis x = 0 || y = 1
function rcvTemp( T1, c_RHS, axis)

T0 = zeros(m_numNos,1)

if c_RHS == 0    # Boundary
    
    if axis == 0
        deltaT = m_nx
        T0[(m_numNos-m_ny):m_numNos] .= deltaT
    elseif axis == 1
        deltaT = m_ny
        T0[1:(m_ny+1):m_numNos] .= deltaT
    end
    T = T1 + T0
    
elseif c_RHS == 1 # Domain
    
    if axis == 0
        flux = 0:m_nx
        T0 = vcat(fill.(flux, (m_ny+1))...)
    elseif axis == 1
        flux = m_ny:-1:0
        T0 = vcat(vcat(fill.(flux, (m_ny+1))'...)...)
    end
    T = T0 - T1
end

return T
end
# -----------------------------------------------------

# -----------------------------------------------------
# Compute Flux - FEM
# Effective property
function femEffective(Tx,Ty)

QXx = 0; QXy  = 0; QYx = 0; QYy = 0 ; V = 0;
for e=1:m_numElem
    c = floor(Int, ((e-1)/m_ny) + 1)
    N1 = e + c;
    N2 = e + c + 1 + m_ny;
    N3 = e + c + m_ny;
    N4 = e + c - 1;
    
    pElemDOFNum = [N1 N2 N3 N4];
    b = m_B[:,:,m_materials[m_matID[e],1],1];

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

# -----------------------------------------------------

# -----------------------------------------------------









# Read File | RAW
# io = open("test6.txt", "r")

# bin_array = read(io)
# result = reshape(reinterpret(UInt8, bin_array), :)

# Read File | NF

m_nx = 3
m_ny = 3

# Initial variables
m_tol = 0
m_numMat = 0
m_numNos = (m_nx+1)*(m_ny+1)
m_numElem = (m_nx)*(m_ny)
m_materials = zeros(UInt64,256,2)

# --------------------
# Lopper over file
m_materials[1,:] = [1 1]
m_materials[256,:] = [2 10]
m_numMat = 2 
m_tol = 10e-6
m_matID = [1;1;1;1;256;1;1;1;1]

# --------------------

# Initial Variables
m_gdlNo = 1
m_nGDL = m_gdlNo*m_nx*m_ny
m_k = zeros(m_gdlNo*4,m_gdlNo*4,m_numMat) 
m_B = zeros(m_gdlNo*2,m_gdlNo*4,m_numMat) 


# Conductivity Matrix Materials
mat = m_materials[:,1]
i = 0
for prop in m_materials[mat.!=0,2]
    global i += 1
    elems_prop = prop
    m_k[:,:,i] , m_B[:,:,i] = Q4ElementConductivity(elems_prop)
end

# Degree of Freedom Map
m_dofMap = dof_modes()
# println(m_dofMap)

# Compute RHS - Boundary or Domain
# c_RHS Boundary = 0 || Domain = 1
#[m_RHS]   % 32 or 64 bits * m_nGDL double  || axis 0 = X
# [m_RHS_2] % 32 or 64 bits * m_nGDL double  || axis 1 = Y
c_RHS = 0
m_RHS   = computeRHS(m_dofMap, c_RHS, 0)
m_RHS_2 = computeRHS(m_dofMap, c_RHS, 1)
# println(m_RHS)
# println(m_RHS_2)


# Preconditioned Conjugate Gradient
# tx 32 or 64 bits * m_nGDL double || ty 64 bits * m_nGDL double
tx = pcg(m_dofMap, m_RHS)
ty = pcg(m_dofMap, m_RHS_2)
# println(tx)
# println(ty)

# Recover Values
# Tx 32 or 64 bits * m_nGDL double || Ty 64 bits * m_nGDL double
Tx = rcvTemp(tx[m_dofMap], c_RHS, 0)
Ty = rcvTemp(ty[m_dofMap], c_RHS, 1)
# println(Tx)
# println(Ty)

# Compute Effective Property
C = femEffective(Tx,Ty)
println(C)