import "boundaryConditions"


type index = {x: i64, y: i64, z: i64}

-- negative z
--      x
--      +  -
-- y +  1  0
--   -  2  3

-- positive z
--      x
--      +  -
-- y +  5  4
--   -  6  7

def getLocalNodeIndex (node :index) :i32 =
  match (node.x,node.y,node.z)
    case (0,1,0) ->  0
    case (1,1,0) ->  1
    case (1,0,0) ->  2
    case (0,0,0) ->  3
    case (0,1,1) ->  4
    case (1,1,1) ->  5
    case (1,0,1) ->  6
    case (0,0,1) ->  7
    case _       -> -1

def getNodeIndices (elementIndex :index) :[8]index =
  let ei = elementIndex
  in [{x=ei.x+0,y=ei.y+1,z=ei.z+0},{x=ei.x+1,y=ei.y+1,z=ei.z+0},
      {x=ei.x+1,y=ei.y+0,z=ei.z+0},{x=ei.x+0,y=ei.y+0,z=ei.z+0},
      {x=ei.x+0,y=ei.y+1,z=ei.z+1},{x=ei.x+1,y=ei.y+1,z=ei.z+1},
      {x=ei.x+1,y=ei.y+0,z=ei.z+1},{x=ei.x+0,y=ei.y+0,z=ei.z+1}]

def getElementIndices (nodeIndex :index) :[8]index =
  let ni = nodeIndex
  in [{x=ni.x-1,y=ni.y+0,z=ni.z-1},{x=ni.x+0,y=ni.y+0,z=ni.z-1},
      {x=ni.x+0,y=ni.y-1,z=ni.z-1},{x=ni.x-1,y=ni.y-1,z=ni.z-1},
      {x=ni.x-1,y=ni.y+0,z=ni.z+0},{x=ni.x+0,y=ni.y+0,z=ni.z+0},
      {x=ni.x+0,y=ni.y-1,z=ni.z+0},{x=ni.x-1,y=ni.y-1,z=ni.z+0}]

def NodeIndicesToDofIndices (nodeIndices :[8]index) =
  let nodes = replicate 3 nodeIndices |> transpose |> flatten_to (24)
  let dofs  = iota 3 |> replicate 8 |> flatten_to (24)
  in (nodes, dofs)

def addIndices (a :index) (b :index) :index =
  {x=a.x+b.x,y=a.y+b.y,z=a.z+b.z}

def indexIsInside (nelx :i64, nely :i64, nelz :i64) (idx :index) :bool =
  (idx.x >= 0 && idx.y >= 0 && idx.z >= 0 && idx.x < nelx && idx.y < nely && idx.z < nelz)

-- loading methods
def loadNodeDisplacement 't [nx][ny][nz] (u :[nx][ny][nz]t) (ni :index) =
  #[unsafe](u[ni.x,ni.y,ni.z])

def readLocalState 't [nx][ny][nz] (u :[nx][ny][nz][3]t) (elementIndex :index) :[24]t =
  let nodeIndices = getNodeIndices elementIndex
  in map (loadNodeDisplacement u) nodeIndices |> flatten_to 24

def getLocalState [nx][ny][nz] (u :[nx][ny][nz][3]f64) (elementIndex :index) :*[24]f64 =
  if (indexIsInside (nx-1,ny-1,nz-1) elementIndex)
  then
    let (nodeIndices, dofIndices) = elementIndex |> getNodeIndices |> NodeIndicesToDofIndices
    let isOnB = map2 (isOnBoundary (nx,ny,nz)) nodeIndices dofIndices
    let ulocal = readLocalState u elementIndex
    in map2 (\b v -> if b then 0 else v) isOnB ulocal
  else
    replicate 24 0

def getDensityUnsafe [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (elementIndex :index) :f64 =
  #[unsafe](f64.f32 x[elementIndex.x,elementIndex.y,elementIndex.z])

def getDensity [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (elementIndex :index) :f64 =
  if indexIsInside (nelx,nely,nelz) elementIndex
    then getDensityUnsafe x elementIndex
    else 0
