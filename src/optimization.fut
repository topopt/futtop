import "solvers"
import "multigrid"
import "densityFilter"
import "complianceSensitivity"
import "utility"

def getVolume [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) :f32 =
  x |> flatten_to (nelx*nely)
    |> flatten_to (nelx*nely*nelz)
    |> reduce (+) 0

def getNewDensity [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (dc :[nelx][nely][nelz]f32) (dv :[nelx][nely][nelz]f32) (lambda :f32) =
  let move = 0.1f32
  in map3_3d (\xx dc dv -> f32.max 0 (f32.max (xx-move) (f32.min 1 (f32.min (xx+move) (xx*f32.sqrt((-dc)/(dv*lambda))))))) x dc dv

def updateDensignOC [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (dc :[nelx][nely][nelz]f32) (dv :[nelx][nely][nelz]f32) (g :f32) =
  let (l1,l2) = loop (l1,l2) = (0f32,1e9f32) while (((l2 - l1) / (l1 + l2)) > 1e-6) do
    let lmid = 0.5 * (l2 + l1)
    let xtest = getNewDensity x dc dv lmid
    let gt   = map3_3d (\xn xo dv -> (xn-xo)*dv) xtest x dv |> flatten_3d |> reduce (+) 0 |> (+g)
    in if (gt > 0) then
      (lmid,l2)
    else
      (l1,lmid)
  let lmid = 0.5 * (l2 + l1)
  let xnew = getNewDensity x dc dv lmid
  let change = map2_3d (\xn xo -> f32.abs(xn-xo)) xnew x |> flatten_3d |> reduce (f32.max) 0
  in (change,xnew)

def performDesignIteration [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (xPhys :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (volfrac :f32) (rmin :f32) =
  let (compliance,dc) = getComplianceSensitivity xPhys u
  let dc = backwardDensityFilter rmin dc
  let dv = replicate nelx (replicate nely (replicate nelz 1f32))
        |> backwardDensityFilter rmin
  let vol = getVolume xPhys
  let volScaled = vol / (f32.i64 (nelx*nely*nelz))
  let g = vol - (volfrac * (f32.i64 (nelx*nely*nelz)))
  let (change,x) = updateDensignOC x dc dv g
  in (compliance,change,volScaled,x)
