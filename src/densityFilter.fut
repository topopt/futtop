type index = {x: i64, y: i64, z: i64}

-- Check if index is valid.
let isInsideDomain nelx nely nelz (idx :index) :bool =
  idx.x >= 0 && idx.y >= 0 && idx.z >= 0 &&
  idx.x < nelx && idx.y < nely && idx.z < nelz

-- Compute weight given radius, element, and neighbour.
let getFilterWeight rmin (ownIdx :index) (neighIdx :index) =
  f32.max 0 (rmin - f32.sqrt( f32.i64 (
    (ownIdx.x-neighIdx.x)*(ownIdx.x-neighIdx.x) +
    (ownIdx.y-neighIdx.y)*(ownIdx.y-neighIdx.y) +
    (ownIdx.z-neighIdx.z)*(ownIdx.z-neighIdx.z))))

-- Sums the output of computeValue in a local neighbourhood with radius boxRadius of element [i,j,k].
let sumOverNeighbourhood boxRadius nelx nely nelz i j k (computeValue :index -> f32) = 
    let boxSize = 2*boxRadius+1
    in tabulate_3d boxSize boxSize boxSize (\ii jj kk ->
        let neighIdx :index = {x=i+ii-boxRadius,y=j+jj-boxRadius,z=k+kk-boxRadius} 
        in
          if isInsideDomain nelx nely nelz neighIdx 
            then computeValue neighIdx 
            else 0)
      |> map (map f32.sum)
      |> map f32.sum
      |> f32.sum

-- given origin and neighbour indices, compute the weighted density contribution.
let getScaledDensity (x :[][][]f32) rmin ownIdx neighIdx = 
      let wgt = getFilterWeight rmin ownIdx neighIdx
      let dens = #[unsafe] x[neighIdx.x,neighIdx.y,neighIdx.z]
      in  (wgt * dens)

-- For all elements, compute the filtered density.
entry forwardDensityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  let boxRadius = i64.max 0 (i64.f32 (f32.ceil (rmin-1)))
  let sumOverThisNeighbourhood = sumOverNeighbourhood boxRadius nelx nely nelz
  in tabulate_3d nelx nely nelz (\i j k->
    let ownIdx :index = {x=i,y=j,z=k}
    let weightSum = sumOverThisNeighbourhood i j k (getFilterWeight rmin ownIdx)
    let scaledDensity = sumOverThisNeighbourhood i j k (getScaledDensity x rmin ownIdx)
    in scaledDensity / weightSum)


entry backwardDensityFilter [nelx][nely][nelz] (rmin :f32) (x :[nelx][nely][nelz]f32) :[nelx][nely][nelz]f32 =
  let boxRadius = i64.max 0 (i64.f32 (f32.ceil (rmin-1)))
  let sumOverThisNeighbourhood = sumOverNeighbourhood boxRadius nelx nely nelz

  let weigthSums =
    tabulate_3d nelx nely nelz (\i j k->
      let ownIdx :index = {x=i,y=j,z=k}
      in sumOverThisNeighbourhood i j k (getFilterWeight rmin ownIdx))

  let getBackwardFilteredDensity ownIdx neighIdx = 
    let weight = getFilterWeight rmin ownIdx neighIdx
    let weightSum = #[unsafe] weigthSums[neighIdx.x,neighIdx.y,neighIdx.z]
    let dens = #[unsafe] x[neighIdx.x,neighIdx.y,neighIdx.z]
    in weight * dens/weightSum

  in tabulate_3d nelx nely nelz (\i j k->
    let ownIdx :index = {x=i,y=j,z=k}
    in sumOverThisNeighbourhood i j k (getBackwardFilteredDensity ownIdx))


-- ==
-- entry: forwardDensityFilter
-- nobench input @../testData/densityFilter.txt output @../testData/densityFilterForwardOut.txt
-- compiled random input { 1.5f32 [64][64][64]f32 }
-- compiled random input { 2.5f32 [64][64][64]f32 }
-- compiled random input { 5.5f32 [64][64][64]f32 }
-- compiled random input { 10.5f32 [64][64][64]f32 }
-- compiled random input { 1.5f32 [128][128][128]f32 }

-- ==
-- entry: backwardDensityFilter
-- nobench input @../testData/densityFilter.txt output @../testData/densityFilterBackwardOut.txt
-- compiled random input { 1.5f32 [64][64][64]f32}
-- compiled random input { 2.5f32 [64][64][64]f32}
-- compiled random input { 5.5f32 [64][64][64]f32}
-- compiled random input { 1.5f32 [128][128][128]f32}
