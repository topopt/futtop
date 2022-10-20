import "keConstants"
import "indexUtilities"
import "boundaryConditions"
import "projection"
import "utility"
import "material"

def getElementContribution [nx][ny][nz][nelx][nely][nelz] (nodeIndex :index, u :[nx][ny][nz][3]f64, x :[nelx][nely][nelz]f32) (elementOffset :index) :[3]f64 =
  let elementIndex = addIndices nodeIndex elementOffset
  let ulocal = getLocalState u elementIndex
  let E = getElementYoungsModule x elementIndex
  let recievingNodeOffset :index = {x=(-elementOffset.x), y=(-elementOffset.y), z=(-elementOffset.z)}
  let li = i64.i32 (getLocalNodeIndex recievingNodeOffset)
  in map (\ke_row -> (f64.sum (map2(\ke u -> ke*E*u) ke_row ulocal))) keslices[li]

def applyStencilOnNode [nelx][nely][nelz][nx][ny][nz] (nodeIndex :index, u :[nx][ny][nz][3]f64, x :[nelx][nely][nelz]f32) :[3]f64 =
   [{x=( 0),y=( 0),z=( 0)}, {x=(-1),y=( 0),z=( 0)},
    {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)},
    {x=( 0),y=( 0),z=(-1)}, {x=(-1),y=( 0),z=(-1)},
    {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)}]
    |> map (\eo -> getElementContribution (nodeIndex, u, x) eo)
    |> transpose
    |> map (\x -> f64.sum x)

-- returns the matrix-vector product K(x)*u
#[noinline]
entry applyStiffnessMatrix [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f64) :[nx][ny][nz][3]f64 =
  let f = tabulate_3d nx ny nz (\i j k -> applyStencilOnNode({x=i,y=j,z=k}, u, x))
  in setBCtoInput u f

-- returns the matrix-vector product K(x)*u for f32
def applyStiffnessMatrixSingle [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32) :[nx][ny][nz][3]f32 =
  map_4d f64.f32 u
  |> applyStiffnessMatrix x
  |> map_4d f32.f64

-- returns the matrix-vector product (P^T*K(x)*P)*u for any subspace
def applyCoarseStiffnessMatrix [nelx][nely][nelz][nxc][nyc][nzc] (l :u8) (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) :[nxc][nyc][nzc][3]f64 =
  let ufine =
    loop u for i < (i16.u8 l) do
      projectToFiner u

  let rfine = applyStiffnessMatrix x ufine

  let rcoarse =
    loop r = rfine for i < (i16.u8 l) do
      projectToCoarser r

  in (rcoarse :> [nxc][nyc][nzc][3]f64)

-- ==
-- entry: applyStiffnessMatrix
-- nobench input @../testData/applyStateOperator1.txt output @../testData/applyStateOperator1Output.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f64 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f64 }
-- compiled random input { [256][256][256]f32 [257][257][257][3]f64 }

entry applyCoarseStiffnessMatrix_test0 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 0 x u

entry applyCoarseStiffnessMatrix_test1 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 1 x u

entry applyCoarseStiffnessMatrix_test2 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 2 x u

entry applyCoarseStiffnessMatrix_test3 [nelx][nely][nelz][nxc][nyc][nzc] (x :[nelx][nely][nelz]f32) (u :[nxc][nyc][nzc][3]f64) =
  applyCoarseStiffnessMatrix 3 x u

-- ==
-- entry: applyCoarseStiffnessMatrix_test0
-- compiled random input { [128][128][128]f32 [129][129][129][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test1
-- compiled random input { [128][128][128]f32 [65][65][65][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test2
-- compiled random input { [128][128][128]f32 [33][33][33][3]f64 }

-- ==
-- entry: applyCoarseStiffnessMatrix_test3
-- compiled random input { [128][128][128]f32 [17][17][17][3]f64 }
