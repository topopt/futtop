import "utility"
import "indexUtilities"
import "assembly"
import "boundaryConditions"
import "material"
import "keConstants"
import "applyStiffnessMatrix"

def omega_const         :f64 = 0.6

def sorLocalMatrix (omega :f64) (f :[3]f64) (S :[3]f64) (M :[3][3]f64) (u :[3]f64) =
  let rx     = #[unsafe] M[0,1]*u[1] + M[0,2]*u[2]
  let ux_new = #[unsafe] (1/M[0,0]) * (f[0]-S[0]-rx)
  let ry     = #[unsafe] M[1,0]*ux_new + M[1,2]*u[2]
  let uy_new = #[unsafe] (1/M[1,1]) * (f[1]-S[1]-ry)
  let rz     = #[unsafe] M[2,0]*ux_new + M[2,1]*uy_new
  let uz_new = #[unsafe] (1/M[2,2]) * (f[2]-S[2]-rz)
  let unew = [ux_new, uy_new, uz_new]
  in #[sequential] map2 (\un uo -> omega*un + (1-omega)*uo) unew u

def sorLocalMatrixBack (omega :f64) (f :[3]f64) (S :[3]f64) (M :[3][3]f64) (u :[3]f64) =
  let rz     = #[unsafe] M[2,0]*u[0] + M[2,1]*u[1]
  let uz_new = #[unsafe] (1/M[2,2]) * (f[2]-S[2]-rz)
  let ry     = #[unsafe] M[1,0]*u[0] + M[1,2]*uz_new
  let uy_new = #[unsafe] (1/M[1,1]) * (f[1]-S[1]-ry)
  let rx     = #[unsafe] M[0,1]*uy_new + M[0,2]*uz_new
  let ux_new = #[unsafe] (1/M[0,0]) * (f[0]-S[0]-rx)
  let unew = [ux_new, uy_new, uz_new]
  in #[sequential] map2 (\un uo -> omega*un + (1-omega)*uo) unew u

def sorForwardAssembled (omega :f64) (mat :[3][81]f64) (uStencil :*[81]f64) (f :[3]f64) :*[3]f64 =
  let uold = #[unsafe] [uStencil[39], uStencil[40], uStencil[41]]
  let uStencil = uStencil with [39] = 0
  let uStencil = uStencil with [40] = 0
  let uStencil = uStencil with [41] = 0
  let S = map (\row -> (f64.sum (map2 (*) uStencil row))) mat
  let M = mat[:3,39:42] :> [3][3]f64
  in sorLocalMatrix omega f S M uold

def sorBackwardAssembled (omega :f64) (mat :[3][81]f64) (uStencil :*[81]f64) (f :[3]f64) :*[3]f64 =
  let uold = #[unsafe] [uStencil[39], uStencil[40], uStencil[41]]
  let uStencil = uStencil with [39] = 0
  let uStencil = uStencil with [40] = 0
  let uStencil = uStencil with [41] = 0
  let S = map (\row -> (f64.sum (map2 (*) uStencil row))) mat
  let M = mat[:3,39:42] :> [3][3]f64
  in sorLocalMatrixBack omega f S M uold

def ssorSweepAssembled [nx][ny][nz] (omega :f64) (mat :[nx][ny][nz][3][81]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let uhalf = tabulate_3d nx ny nz (\i j k ->
    let uloc = getLocalU u {x=i,y=j,z=k}
    in #[unsafe] sorForwardAssembled omega mat[i,j,k] uloc f[i,j,k])
  let unew = tabulate_3d nx ny nz (\i j k ->
    let uloc = getLocalU uhalf {x=i,y=j,z=k}
    in #[unsafe] sorBackwardAssembled omega mat[i,j,k] uloc f[i,j,k])
  in unew

entry sorAssembled [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let number_sweeps :i32 = 2
  in (iterate number_sweeps (ssorSweepAssembled omega_const mat f)) u

entry symGS [nx][ny][nz] (mat :[nx][ny][nz][3][81]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let number_sweeps :i32 = 1
  in (iterate number_sweeps (ssorSweepAssembled 1.0 mat f)) u


-- let generateNodeOffsetWithoutCenter (eo :index) =
--   let no = getNodeIndices eo
--   let no = filter (\i -> i.x != 0 || i.y != 0 || i.z != 0) no
--   let no = no :> [7]index
--   in zip (replicate 7 eo) no
--
-- let allOffsetPairs =
--   [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
--    {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
--    {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
--   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}] |> map generateNodeOffsetWithoutCenter |> flatten_to (7*8)

def allOffsetPairs = [
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = 0i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = -1i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = -1i64, y = 1i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = 0i64, y = 1i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64},{x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = 0i64, z = 0i64}, {x = -1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 1i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 0i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 1i64, y = 0i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64}, {x = 1i64, y = -1i64, z = 1i64}),
({x = 0i64, y = -1i64, z = 0i64},{x = 0i64, y = -1i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = 0i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = 0i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = 0i64, y = -1i64, z = 1i64}),
({x = -1i64, y = -1i64, z = 0i64}, {x = -1i64, y = -1i64, z = 1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 1i64, y = 1i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 1i64, y = 0i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 0i64, y = 0i64, z = -1i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = -1i64},{x = 1i64, y = 1i64, z = 0i64}),
({x = 0i64, y = 0i64, z = -1i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64},{x = -1i64, y = 1i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 0i64, z = -1i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = 0i64, y = 1i64, z = 0i64}),
({x = -1i64, y = 0i64, z = -1i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = 0i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = -1i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = -1i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 1i64, y = 0i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64},{x = 1i64, y = -1i64, z = 0i64}),
({x = 0i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = 0i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = 0i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = -1i64, z = -1i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = 0i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = 0i64, y = -1i64, z = 0i64}),
({x = -1i64, y = -1i64, z = -1i64}, {x = -1i64, y = -1i64, z = 0i64})]

def getSendingNode (elementOffset :index, nodeOffset :index) :(i32, i32) =
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let sendingNodeOffset :index = {x=nodeOffset.x-elementOffset.x, y=nodeOffset.y-elementOffset.y, z=nodeOffset.z-elementOffset.z}
  in (getLocalNodeIndex(recievingNodeOffset), getLocalNodeIndex(sendingNodeOffset))

def getInputVector [nx][ny][nz] (nodeIndex :index, nodeOffset :index, u :[nx][ny][nz][3]f64) :[3]f64 =
  let loadIndex :index = {x=nodeIndex.x+nodeOffset.x,y=nodeIndex.y+nodeOffset.y,z=nodeIndex.z+nodeOffset.z} in
  tabulate 3 (\d ->
    if ((indexIsInside (nx,ny,nz) loadIndex) && !(isOnBoundary (nx,ny,nz) loadIndex d)) then
      #[unsafe] (u[loadIndex.x,loadIndex.y,loadIndex.z,d])
    else
      0)

def getInputVectorSingle [nx][ny][nz] (nodeIndex :index, nodeOffset :index, u :[nx][ny][nz][3]f32) :[3]f32 =
  let loadIndex :index = {x=nodeIndex.x+nodeOffset.x,y=nodeIndex.y+nodeOffset.y,z=nodeIndex.z+nodeOffset.z} in
  tabulate 3 (\d ->
    if ((indexIsInside (nx,ny,nz) loadIndex) && !(isOnBoundary (nx,ny,nz) loadIndex d)) then
      #[unsafe] (u[loadIndex.x,loadIndex.y,loadIndex.z,d])
    else
      0)

def multiplyScaledLocalMatrix(m :localMatrix, a :[3]f64, s :f64) :[3]f64 =
  [(s*(m.xx*a[0]+m.xy*a[1]+m.xz*a[2])),
   (s*(m.yx*a[0]+m.yy*a[1]+m.yz*a[2])),
   (s*(m.zx*a[0]+m.zy*a[1]+m.zz*a[2]))]

def scaleLocalMatrix(m :localMatrix) (s :f64) :localMatrix =
 {xx=s*m.xx,xy=s*m.xy,xz=s*m.xz,
  yx=s*m.yx,yy=s*m.yy,yz=s*m.yz,
  zx=s*m.zx,zy=s*m.zy,zz=s*m.zz}

def addLocalMatrix(a :localMatrix) (b :localMatrix) :localMatrix =
 {xx=a.xx+b.xx,xy=a.xy+b.xy,xz=a.xz+b.xz,
  yx=a.yx+b.yx,yy=a.yy+b.yy,yz=a.yz+b.yz,
  zx=a.zx+b.zx,zy=a.zy+b.zy,zz=a.zz+b.zz}

def getLocalMatrix (elementOffset :index, nodeOffset :index) :localMatrix =
  let (recieve,send) = getSendingNode(elementOffset, nodeOffset)
  in getke_l0 (recieve,send)

def getSContribution [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (nodeIndex :index) (elementOffset :index, nodeOffset :index) :[3]f64 =
  let elementIndex = addIndices nodeIndex elementOffset
  let elementScale = getElementYoungsModule x elementIndex
  let localMatrix  = getLocalMatrix(elementOffset,nodeOffset)
  let ulocal       = getInputVectorSingle(nodeIndex,nodeOffset,u)
  let ulocal       = #[sequential] map f64.f32 ulocal
  in multiplyScaledLocalMatrix(localMatrix, ulocal, elementScale)

def build_S [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (nodeIndex :index) :[3]f64 =
  copy allOffsetPairs
  |> (#[sequential] map (getSContribution x u nodeIndex))
  |> transpose
  |> (#[sequential] map f64.sum)

def getOwnLocalMatrix (elementOffset :index) =
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let li = getLocalNodeIndex(recievingNodeOffset)
  in getke_l0 (li,li)

def getMContribution [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (nodeIndex :index) (elementOffset :index) =
  let elementIndex = addIndices nodeIndex elementOffset
  let elementScale = getElementYoungsModule x elementIndex
  let localMatrix  = getOwnLocalMatrix elementOffset
  in scaleLocalMatrix localMatrix elementScale

def build_M [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (nodeIndex :index) =
  [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
   {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
   {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}]
   |> (#[sequential] map (getMContribution x nodeIndex))
   |> (#[sequential] reduce addLocalMatrix {xx=0,xy=0,xz=0,yx=0,yy=0,yz=0,zx=0,zy=0,zz=0})

def sorLocal (f :[3]f64) (S :[3]f64) (M :localMatrix) (u :[3]f64) =
  let rx     = #[unsafe] M.xy*u[1] + M.xz*u[2]
  let ux_new = #[unsafe] (1/M.xx) * (f[0]-S[0]-rx)
  let ry     = #[unsafe] M.yx*ux_new + M.yz*u[2]
  let uy_new = #[unsafe] (1/M.yy) * (f[1]-S[1]-ry)
  let rz     = #[unsafe] M.zx*ux_new + M.zy*uy_new
  let uz_new = #[unsafe] (1/M.zz) * (f[2]-S[2]-rz)
  in [ux_new, uy_new, uz_new]

def sorLocalBack (f :[3]f64) (S :[3]f64) (M :localMatrix) (u :[3]f64) =
  let rz     = #[unsafe] M.zx*u[0] + M.zy*u[1]
  let uz_new = #[unsafe] (1/M.zz) * (f[2]-S[2]-rz)
  let ry     = #[unsafe] M.yx*u[0] + M.yz*uz_new
  let uy_new = #[unsafe] (1/M.yy) * (f[1]-S[1]-ry)
  let rx     = #[unsafe] M.xy*uy_new + M.xz*uz_new
  let ux_new = #[unsafe] (1/M.xx) * (f[0]-S[0]-rx)
  in [ux_new, uy_new, uz_new]

def sorNodeForward  [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (f :[3]f32) (nodeIndex :index) :[3]f32 =
 -- extract value of own node, and zero
 let ni = nodeIndex

 let ux_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,0]
 let uy_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,1]
 let uz_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,2]
 let uold = [ux_old, uy_old, uz_old]

 let f = #[sequential] map f64.f32 f
 let S = build_S x u ni
 let M = build_M x ni

 let utmp = sorLocal f S M uold
 let usmoothed = #[sequential] map2 (\un uo -> omega_const*un + (1-omega_const)*uo) utmp uold
 in #[sequential] map f32.f64 usmoothed

def sorNodeBackward [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) (f :[3]f32) (nodeIndex :index) :[3]f32 =
 -- extract value of own node, and zero
 let ni = nodeIndex

 let ux_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,0]
 let uy_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,1]
 let uz_old = #[unsafe] f64.f32 u[ni.x,ni.y,ni.z,2]
 let uold = [ux_old, uy_old, uz_old]

 let f = #[sequential] map f64.f32 f
 let S = build_S x u ni
 let M = build_M x ni

 let utmp = sorLocalBack f S M uold
 let usmoothed = #[sequential] map2 (\un uo -> omega_const*un + (1-omega_const)*uo) utmp uold
 in #[sequential] map f32.f64 usmoothed

def ssorSweep [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
  let uhalf = tabulate_3d nx ny nz (\i j k -> sorNodeForward x u f[i,j,k] {x=i,y=j,z=k})
  let uhalf = setBCtoZero 0 uhalf
  let unew = tabulate_3d nx ny nz (\i j k -> sorNodeBackward x uhalf f[i,j,k] {x=i,y=j,z=k})
  in setBCtoZero 0 unew

entry sorMatrixFree [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (f :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
  let number_sweeps :i32 = 1
  in (iterate number_sweeps (ssorSweep x f)) u

-- entry jacobiReference [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (b :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
--  let invD = assembleInverseDiagonal 0 x |> map_4d f32.f64 :> [nx][ny][nz][3]f32
--  in jacobiSmootherSingle x b invD u

-- ==
-- entry: sorMatrixFree
-- nobench input @../testData/sor1.txt auto output
-- nobench input @../testData/sor1.txt output @../testData/sor1Out.txt
-- nobench input @../testData/sor2.txt output @../testData/sor2Out.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f32 [65][65][65][3]f32 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f32 [129][129][129][3]f32 }
-- compiled random input { [256][256][256]f32 [257][257][257][3]f32 [257][257][257][3]f32 }

def sorCoarseForward [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let Ku = applyCoarseStiffnessMatrix l x u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d (sorLocalMatrix omega_const) f S bdiag u
  in setBCtoZero 0 unew

def sorCoarseBackward [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let Ku = applyCoarseStiffnessMatrix l x u
  let Kd = map2_3d vecmul_f64 bdiag u
  let S  = map2_4d (-) Ku Kd
  let unew = map4_3d (sorLocalMatrixBack omega_const) f S bdiag u
  in setBCtoZero 0 unew

def sorCoarseSweep [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) = 
  let uhalf = sorCoarseForward l x bdiag f u
  in sorCoarseBackward l x bdiag f uhalf

def sorCoarse [nelx][nely][nelz][nx][ny][nz] (l :u8) (x :[nelx][nely][nelz]f32) (bdiag :[nx][ny][nz][3][3]f64) (f :[nx][ny][nz][3]f64) (u :[nx][ny][nz][3]f64) =
  let sweep v = sorCoarseSweep l x bdiag f v
  let number_sweeps :i32 = 1
  in (iterate number_sweeps sweep) u


