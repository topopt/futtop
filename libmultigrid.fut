import "src/solvers"
import "src/multigrid"
import "src/densityFilter"
import "src/optimization"

-- Applies the density filter with radius rmin to the structured grid variables given in x.
entry densityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  forwardDensityFilter rmin x

-- Assemble the multigrid data required to solve the linear system of equations for a given density field x.
entry assembleMultigridData [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :multigridData =
  generateMultigridData x

-- Solve the system of eqations for u K(xPhys) u = f, using uold as initial guess.
entry solveSystem [nelx][nely][nelz][nx][ny][nz] (xPhys :[nelx][nely][nelz]f32) (mgData :multigridData) (uold :[nx][ny][nz][3]f32) (f :[nx][ny][nz][3]f32) =
  cgSolveMG xPhys mgData f uold

-- Update design variables x given displacement field and volume fraction.
entry designIteration [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (xPhys :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (volfrac :f32) (rmin :f32) =
  performDesignIteration x xPhys u volfrac rmin

-- ==
-- entry: assembleMultigridData
-- compiled random input { [40][20][20]f32 }
-- compiled random input { [80][40][40]f32 }
-- compiled random input { [160][80][80]f32 }

-- ==
-- entry: densityFilter
-- compiled random input { 2.5f32 [40][20][20]f32 }
-- compiled random input { 2.5f32 [80][40][40]f32 }
-- compiled random input { 2.5f32 [160][80][80]f32 }
-- compiled random input { 1.5f32 [160][80][80]f32 }
-- compiled random input { 5.5f32 [160][80][80]f32 }

-- ==
-- entry: designIteration
-- compiled random input { [80][40][40]f32    [80][40][40]f32    [81][41][41][3]f32     0.12f32 2.5f32 }
-- compiled random input { [160][80][80]f32   [160][80][80]f32   [161][81][81][3]f32    0.12f32 2.5f32 }
-- compiled random input { [320][160][160]f32 [320][160][160]f32 [321][161][161][3]f32  0.12f32 2.5f32 }
