from firedrake import *
import numpy as np

def check_nullspace(residual_form, solution_variable, null_space_candidate):
    """
    Symbolically checks if a vector `null_space_candidate` lies in the
    nullspace of the Jacobian of a given residual form.
    """
    du = TrialFunction(solution_variable.function_space())
    J_form = derivative(residual_form, solution_variable, du)
    J_z = replace(J_form, {du: null_space_candidate})
    
    # Assemble the 1-form J*z, which results in a Cofunction
    result_vector = assemble(J_z)

    # A Cofunction is not a UFL object. To get its norm,
    # we access the underlying PETSc vector and use its norm method.
    norm_Jz = result_vector.vector().norm()

    print(f"Norm of (Jacobian * null_space_candidate) is: {norm_Jz}")
    return np.isclose(norm_Jz, 0.0)

# --- Main script for Stokes ---

# 1. Set up the problem
mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# The 'state vector' is a mixed function for velocity and pressure
w = Function(W, name="StateVector")
(u, p) = split(w)

# Test functions
(v, q) = TestFunctions(W)

# A parameter for viscosity, just to make it more realistic
nu = Constant(1.0)

# 2. Write the global residual for Stokes
F = (
    dot(dot(grad(u), u), v) * dx
    + nu * inner(sym(grad(u)), sym(grad(v))) * dx
    - p * div(v) * dx
    + q * div(u) * dx
)

# 3. Define a candidate for the null space: (u=0, p=1)
null_space_candidate = Function(W, name="NullSpaceCandidate")
u_ns_cand, p_ns_cand = null_space_candidate.subfunctions
p_ns_cand.assign(1.0) 

# 4. Check for the nullspace
print("Checking the nonlinear Stokes problem:")
is_singular = check_nullspace(F, w, null_space_candidate)

if is_singular:
    print("The Stokes Jacobian is singular for a constant pressure vector. ✅")
    print("The system is not closed and needs a pressure constraint.")
else:
    print("The Stokes Jacobian is not singular for the given vector.")