import "applyStiffnessMatrix"
import "assembly"
import "projection"
import "utility"
import "sor"

type~ mgL3Data = ([][][][3]f64, [][][][3][81]f64)
type~ mgL2Data = [][][][3][81]f64
type~ mgL1Data = {}
type~ mgL0Data = {}
type~ multigridData = (mgL0Data, mgL1Data, mgL2Data, mgL3Data)

def generateMultigridData [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :multigridData =
  -- let d1 = assembleBlockDiagonal 1 x
  let m2 = assembleStiffnessMatrix 2 x
  let m3 = assembleStiffnessMatrix 3 x
  let d3 = extractInverseDiagonal m3
  in ({}, {}, m2, (d3,m3))

def jacobiSmoother [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (invD :[nx][ny][nz][3]f64) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let nsweeps = 4
  let omega   = 0.6
  let smooth u = applyCoarseStiffnessMatrix l x u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u

def jacobiSmootherFine [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (invD :[nx][ny][nz][3]f64) (b :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  let nsweeps = 4
  let omega   = 0.6
  let b = map_4d f64.f32 b
  let u = map_4d f64.f32 u
  let smooth u = applyStiffnessMatrix x u
    |> map4_4d (\uu bb dd tt -> uu - omega * dd * (tt - bb)) u b invD
  in (iterate nsweeps smooth) u |> map_4d f32.f64

def jacobiBlockSmoother [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (b :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let nsweeps = 6
  let omega   = 0.6
  let smooth u = applyCoarseStiffnessMatrix l x u
    |> map2_4d (\bb rr -> rr - bb) b
    |> applyBlockDiagonal bdiag
    |> map2_4d (\uold unew -> uold - omega * unew) u
  in (iterate nsweeps smooth) u

def innerProduct [n][m][l][k] (a :[n][m][l][k]f64) (b :[n][m][l][k]f64) :f64 =
  map2_4d (*) a b
  |> map_3d f64.sum
  |> map_2d f64.sum
  |> map    f64.sum
  |>        f64.sum

def norm [n][m][l][k] (a :[n][m][l][k]f64) :f64 =
  map_4d (\x -> x*x) a
  |> map_3d f64.sum
  |> map_2d f64.sum
  |> map    f64.sum
  |>        f64.sum
  |> f64.sqrt

def cgSolveJacSubspace [nx][ny][nz] (data :multigridData) (b :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let maxIt   = 800
  let invD    = data.3.0 :> [nx][ny][nz][3]f64
  let matrix  = data.3.1 :> [nx][ny][nz][3][81]f64
  let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))

  let inner_iteration (uold, rold, pold, rhoold) = 
    let z      = map2_4d (*) invD rold -- jacobi preconditioner
    let rho    = innerProduct rold z
    let beta   = rho / rhoold
    let p      = map2_4d (\pp zz -> beta * pp + zz) pold z
    let q      = applyAssembledStiffnessMatrix matrix p
    let alpha  = rho / (innerProduct p q)
    let u      = map2_4d (\uu pp -> uu + alpha * pp) uold p
    let r      = map2_4d (\rr qq -> rr - alpha * qq) rold q
    in (u, r, p, rho)

  let  (u, _, _, _) = (iterate maxIt inner_iteration) (zero_4d, b, zero_4d, 1f64)
  in u

def vcycle_l2 [nx][ny][nz] (data :multigridData) (f :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let matrix = data.2 :> [nx][ny][nz][3][81]f64
  let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
  let z = sorAssembled matrix f zero_4d
  let v = z
    |> applyAssembledStiffnessMatrix matrix
    |> map2_4d (\ff dd -> ff - dd) f
    |> projectToCoarser
    |> cgSolveJacSubspace data
    |> projectToFiner
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
  in sorAssembled matrix f z

-- def vcycle_l1 [nelx][nely][nelz][nx][ny][nz] (data :multigridData) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
--   let diag = data.1 :> [nx][ny][nz][3][3]f64
--   let zero_4d = replicate nx (replicate ny (replicate nz [0f64,0f64,0f64]))
--   let z = sorCoarse 1 x diag f zero_4d
--   let v = z
--     |> applyCoarseStiffnessMatrix 1 x
--     |> map2_4d (\ff dd -> ff - dd) f
--     |> projectToCoarser
--     |> vcycle_l2 data
--     |> projectToFiner 
--   let z = map2_4d (+) z (v :> [nx][ny][nz][3]f64)
--   in sorCoarse 1 x diag f z

def vcycle_l0 [nelx][nely][nelz][nx][ny][nz] (data :multigridData) (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  let zero_4d = replicate nx (replicate ny (replicate nz [0,0,0f32]))
  let z = sorMatrixFree x f zero_4d
  let v = z
    |> applyStiffnessMatrixSingle x
    |> map2_4d (\ff dd -> ff - dd) f -- compute r
    |> projectToCoarserSingle
    |> projectToCoarserSingle
    |> map_4d f64.f32
    |> vcycle_l2 data -- compute v
    |> map_4d f32.f64
    |> projectToFinerSingle
    |> projectToFinerSingle
  let z = map2_4d (+) z (v :> [nx][ny][nz][3]f32)
  in sorMatrixFree x f z

