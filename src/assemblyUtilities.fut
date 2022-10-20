import "assemblyWeights"
import "indexUtilities"
import "boundaryConditions"
import "utility"

type localState   = [24]f64
type nodalWeights = [8]f64

def elementOffsets :[8]index =
  [{x=(-1),y=( 0),z=(-1)}, {x=( 0),y=( 0),z=(-1)},
   {x=( 0),y=(-1),z=(-1)}, {x=(-1),y=(-1),z=(-1)},
   {x=(-1),y=( 0),z=( 0)}, {x=( 0),y=( 0),z=( 0)},
   {x=( 0),y=(-1),z=( 0)}, {x=(-1),y=(-1),z=( 0)}]

-- generates initial weights with zeros and 1 on the local index
def getInitialWeights (nodeOffset :index) :nodalWeights =
 let li = getLocalNodeIndex nodeOffset
 in (replicate 8 0) with [li] = 1

-- zeroes out indices on boundary
def applyBoundaryConditionsToWeights (elementIndex :index) (w :nodalWeights) (nx :i64, ny :i64, nz :i64) (dofNumber :i64) :nodalWeights =
 getNodeIndices elementIndex
   |> (#[sequential] map (\x -> isOnBoundary (nx,ny,nz) x dofNumber))
   |> (#[sequential] map2 (\v b -> if b then 0 else v) w)

-- zeroes out indices not on boundary
def applyBoundaryConditionsToWeightsInverse (elementIndex :index) (w :nodalWeights) (nx :i64, ny :i64, nz :i64) (dofNumber :i64) :nodalWeights =
 getNodeIndices elementIndex
   |> (#[sequential] map (\x -> isOnBoundary (nx,ny,nz) x dofNumber && (indexIsInside (nx,ny,nz) x)))
   |> (#[sequential] map2 (\v b -> if b then (v/8) else 0) w)

def generateLoad (o :i64) (w :nodalWeights) :*localState =
 #[sequential] scatter (replicate 24 0) [0+o,3+o,6+o,9+o,12+o,15+o,18+o,21+o] w

def prolongateCellIndices (cellIndex :index) =
  getNodeIndices {x=2*cellIndex.x, y=2*cellIndex.y, z=2*cellIndex.z}

-- prolongs a set of eight nodal weights and the cell index
def prolongateCellValues (w :nodalWeights) :[8]nodalWeights =
  #[sequential] map (\pMat -> vecmul_f64 pMat w) prolongationWeights

-- combines eight values to one. Is the transpose operation of prolongateCellValues
def restrictCell (vals :[8][24]f64) :[24]f64 =
  let indX   = [0,3,6, 9,12,15,18,21]
  let indY   = [1,4,7,10,13,16,19,22]
  let indZ   = [2,5,8,11,14,17,20,23]
  let indAll = [0,3,6, 9,12,15,18,21,1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23]

  let valX = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indX) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  let valY = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indY) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  let valZ = #[sequential] map (\v -> map (\i -> #[unsafe](v[i])) indZ) vals
    |> map2 (\pMat a -> vecmul_f64 (transpose pMat) a) prolongationWeights
    |> transpose
    |> map (reduce (+) 0)

  in scatter (replicate 24 0) indAll ([valX,valY,valZ] |> flatten_to 24)

-- computes the contribution of a coarse grid cell to a
def getCoarseCellContribution (getFineValue: index -> nodalWeights -> localState) (l :u8) (cellIndex :index) (w :nodalWeights) =
  -- this is implemented as a case-match to allow the compiler to reason about
  -- array sizes, and realize that we will always restrict to an array of 24
  -- values.
  match l
    case 0 ->
      getFineValue cellIndex w
    case 1 ->
      let fineCellWeights = prolongateCellValues w
      let fineCellIndices = prolongateCellIndices cellIndex
      let fineValues = map2 getFineValue fineCellIndices fineCellWeights
      in restrictCell fineValues
    case 2 ->
      let fineCellWeights = prolongateCellValues w
        |> map prolongateCellValues
      let fineCellIndices = prolongateCellIndices cellIndex
        |> map prolongateCellIndices
      let fineValues = map2_2d getFineValue fineCellIndices fineCellWeights
      in fineValues
        |> map restrictCell
        |>     restrictCell
    case 3 ->
      let fineCellWeights = prolongateCellValues w
        |> map    prolongateCellValues
        |> map_2d prolongateCellValues
      let fineCellIndices = prolongateCellIndices cellIndex
        |> map    prolongateCellIndices
        |> map_2d prolongateCellIndices
      let fineValues = map2_3d getFineValue fineCellIndices fineCellWeights
      in fineValues
        |> map_2d restrictCell
        |> map    restrictCell
        |>        restrictCell
    case 4 ->
      let fineCellWeights = prolongateCellValues w
        |> map    prolongateCellValues
        |> map_2d prolongateCellValues
        |> map_3d prolongateCellValues
      let fineCellIndices = prolongateCellIndices cellIndex
        |> map    prolongateCellIndices
        |> map_2d prolongateCellIndices
        |> map_3d prolongateCellIndices
      let fineValues = map2_4d getFineValue fineCellIndices fineCellWeights
      in fineValues
        |> map_3d restrictCell
        |> map_2d restrictCell
        |> map    restrictCell
        |>        restrictCell
    case _ ->
      replicate 24 0 -- not implemented

-- takes an offset with values [-1..1], and maps each compination to the three local indices
-- let offsetToAssembledOffset (idx :index) =
--   let ind = (idx.x+1)*9 + (idx.y+1)*3 + (idx.z+1)
--   in [3*ind+0,3*ind+1,3*ind+2]
--
-- let elementAssembledOffsets :[192]i64 = elementOffsets
--                                         |> map getNodeIndices
--                                         |> flatten_to 64
--                                         |> map offsetToAssembledOffset
--                                         |> flatten_to 192

def elementAssembledOffsets :[192]i64 =
  [18i64, 19i64, 20i64, 45i64, 46i64, 47i64, 36i64, 37i64, 38i64, 9i64, 10i64,
  11i64, 21i64, 22i64, 23i64, 48i64, 49i64, 50i64, 39i64, 40i64, 41i64, 12i64,
  13i64, 14i64, 45i64, 46i64, 47i64, 72i64, 73i64, 74i64, 63i64, 64i64, 65i64,
  36i64, 37i64, 38i64, 48i64, 49i64, 50i64, 75i64, 76i64, 77i64, 66i64, 67i64,
  68i64, 39i64, 40i64, 41i64, 36i64, 37i64, 38i64, 63i64, 64i64, 65i64, 54i64,
  55i64, 56i64, 27i64, 28i64, 29i64, 39i64, 40i64, 41i64, 66i64, 67i64, 68i64,
  57i64, 58i64, 59i64, 30i64, 31i64, 32i64, 9i64, 10i64, 11i64, 36i64, 37i64,
  38i64, 27i64, 28i64, 29i64, 0i64, 1i64, 2i64, 12i64, 13i64, 14i64, 39i64,
  40i64, 41i64, 30i64, 31i64, 32i64, 3i64, 4i64, 5i64, 21i64, 22i64, 23i64,
  48i64, 49i64, 50i64, 39i64, 40i64, 41i64, 12i64, 13i64, 14i64, 24i64, 25i64,
  26i64, 51i64, 52i64, 53i64, 42i64, 43i64, 44i64, 15i64, 16i64, 17i64, 48i64,
  49i64, 50i64, 75i64, 76i64, 77i64, 66i64, 67i64, 68i64, 39i64, 40i64, 41i64,
  51i64, 52i64, 53i64, 78i64, 79i64, 80i64, 69i64, 70i64, 71i64, 42i64, 43i64,
  44i64, 39i64, 40i64, 41i64, 66i64, 67i64, 68i64, 57i64, 58i64, 59i64, 30i64,
  31i64, 32i64, 42i64, 43i64, 44i64, 69i64, 70i64, 71i64, 60i64, 61i64, 62i64,
  33i64, 34i64, 35i64, 12i64, 13i64, 14i64, 39i64, 40i64, 41i64, 30i64, 31i64,
  32i64, 3i64, 4i64, 5i64, 15i64, 16i64, 17i64, 42i64, 43i64, 44i64, 33i64,
  34i64, 35i64, 6i64, 7i64, 8i64]
