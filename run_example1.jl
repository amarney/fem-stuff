include("finite_element_package_1d.jl")
## ----------------------------------------------------------------------------
# Consider the following BVP with mixed b.c.:
## ----------------------------------------------------------------------------
# BVP: Find u such that
#
#   -(au')' + cu = f,  x \in (0, 1)
#   a(0)u'(0) = 2,    u(1) = exp(1)
#
#   where the known functions are
#
#   a = 1 + cos(x),    c = 1/(x + 1)
#   f = (sin(x) - cos(x) - 1 + 1/(x + 1))e^x
## ----------------------------------------------------------------------------
# The solution is u = e^x.
# From this BVP, we can derive the following weak problem:
## ----------------------------------------------------------------------------
# Weak Problem: Find u \in S such that
#
#   \int_0^1 au'v' dx + \int_0^1 cuv dx = \int_0^1 vf dx - 2v(1)
#
#   for all v \in T where
#
#   T = H^1 (0, 1) with v(0) = 0  
#   S = H^1 (0, 1) with u(1) = exp(1)
## ----------------------------------------------------------------------------
# We call T the "test function space" and S the "trial function set".
# From this weak problem, we can derive the finite element problem:
## ----------------------------------------------------------------------------
# FE problem: Find a vector uGlobal such that
#
#   (A_{1,1} + C_{0,0})*uGlobal = f_{0} + b_{eA} + b_{eC} + b_{n}
#
#   where the known matrices and vectors are
#
#   A_{1,1} is the matrix associated with the bilinear form: \int_0^1 au'v' dx
#   C_{0,0} is the matrix associated with the bilinear form: \int_0^1 cuv dx
#   f_{0} is the vector associated with the linear form: \int_0^1 vf dx
#   b_{eA} is the essential boundary condition vector extracted from A_{1,1}
#   b_{eC} is the essential boundary condition vector extracted from C_{0,0}
#   b_{n} is the vector associated with the boundary term: -2v(1)
## ----------------------------------------------------------------------------

## Make a data structure for the mesh
domain = [0.0, 1.0]  # x \in (0, 1)
nEle = 20            # 20 elements mesh
mesh = mesh1d(domain, nEle)

## Make a data structure to access the finite element shape functions
degree = 2  # second order FE space
fem = fem1d(mesh, degree)

## note:
# fem has 4 fields which correspond exactly to the arrays generated by
# Jeff Borggaard's onedMesh function

## Use Galerkin approximation: same finite element space for u and v
ufem = fem
vfem = fem

## Make additional arrays to aid in solving the finite element problem
indexBC = [-1, 0]  # the left endpoint is natural -1, the right is essential 0
dof, nt = boundaryclassify1d(fem, domain, indexBC)

## note:
# if ufem != vfem, then create dof_u, nt_u, seperate from dof_v, nt_v

## Assemble the f_{0} vector
fName(x) = exp(x).*(sin(x) - cos(x) - 1 + 1./(x + 1))  # fName from linear form
derivOrder = 0  # no derivatives on v
rule = 5        # 5 point quadrature rule to be safe
fh_0 = linformglob1d(fName, mesh, vfem, derivOrder, rule)
f_0 = fh_0[dof]

## Assemble the A_{1,1} matrix
aName(x) = 1 + cos(x)  # aName from bilinear form
derivOrders = [1, 1]   # 1 derivative on u, 1 derivative on v
rule = 5               # 5 point quadrature rule to be safe
Ah_11 = biformglob1d(aName, mesh, ufem, vfem, derivOrders, rule)
A_11 = Ah_11[dof, dof]

## Assemble the C_{0,0} matrix
cName(x) = 1./(x + 1)  # cName from bilinear form
derivOrders = [0, 0]   # no derivatives on u or v
rule = 5               # 5 point quadrature rule to be safe
Ch_00 = biformglob1d(cName, mesh, ufem, vfem, derivOrders, rule)
C_00 = Ch_00[dof, dof]

## Assemble the b_{eA} vector
essName(x) = exp(1)      # essential b.c. enforced on solution set S
uEss = zeros(length(fem.x), 1)  # prepare a vector for essential b.c.
indEss = find(nt .== 0)          # find the indices for essential control pnts
uEss[indEss] = essName(fem.x[indEss])
temp = -Ah_11*uEss   # negative because we subtract it over to the rhs
b_eA = temp[dof]

## Assemble the b_{eC} vector
temp = -Ch_00*uEss   # negative because we subtract it over to the rhs
b_eC = temp[dof]

## Assemble the b_{n} vector
natName(x) = -2*ones(size(x))   # natural b.c. vector from -2*v(1)
uNat = zeros(length(fem.x), 1)  # prepare a vector for natural b.c.
indNat = find(nt .== -1)         # find the indices for the natural control pnts
uNat[indNat] = natName(fem.x[indNat])
b_n = uNat[dof]
# b_n[1] = -b_n[1] # for treating natural boundary condition on the left (?)

## Solve the system of equations and assemble uGlobal
uhGlobal = (A_11 + C_00)\(f_0 + b_eA + b_eC + b_n)
uGlobal = uEss
uGlobal[dof] = uhGlobal

## note: uGlobal contains the coefficients for our uFE basis functions
## ----------------------------------------------------------------------------
## Post-processing section:

## Evaluate our finite element solution at at the point 0.5
tx = 0.5
utx = Float64
for k = 1:size(mesh.eleConn, 1)  # loop over all the elements
  element = mesh.x[mesh.eleConn[k, :]]
  if ((tx - element[1])*(tx - element[2]) <= 0)  # if x is in this element
    uLocal = uGlobal[fem.feConn[k, :]]           # extract the local coefficients
    derivOrder = 0
    utx = evalfe1d(tx, uLocal, element, degree, derivOrder)
    break
  end
end

# for example, compare utx with exp(0.5)

## Compute the global error measured by the L^2 norm
uName(x) = exp(x)  # true solution
rule = 5
derivOrder = 0
u0Error = errorglob1d(uName, uGlobal, mesh, fem, derivOrder, rule)

## Compute the global error measured by the H^1 norm
udName(x) = exp(x)  # derivative of true solution
derivOrder = 1
u1Error = errorglob1d(udName, uGlobal, mesh, fem, derivOrder, rule)
