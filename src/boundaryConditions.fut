type index       = {x: i64, y: i64, z: i64}

-- cantilever - classic
-- def isOnBoundary (nx :i64, ny :i64, nz :i64) (nodeIndex :index) (dofnumber :i64) :bool = (nodeIndex.x == 0)

-- cantilever - sides
def isOnBoundary (nx :i64, ny :i64, nz :i64) (nodeIndex :index) (dofnumber :i64) :bool = 
  (nodeIndex.x == 0) && ((nodeIndex.y < (ny/4)) || (nodeIndex.y > ny-1-(ny/4)))

-- MBB beam
-- def isOnBoundary (nx :i64, ny :i64, nz :i64) (nodeIndex :index) (dofnumber :i64) :bool = 
--   (nodeIndex.x == 0 && dofnumber == 0) || ((nodeIndex.x > (nx-1-(nx/10))) && nodeIndex.z == 0 && dofnumber != 0)

-- snapthrough
-- def isOnBoundary (nx :i64, ny :i64, nz :i64) (nodeIndex :index) (dofnumber :i64) :bool = 
--   (nodeIndex.x == 0 && dofnumber == 0) || (nodeIndex.x == (nx-1))

def setBCtoZero [nx][ny][nz][ndof] 'a (zero :a) (v :[nx][ny][nz][ndof]a) :[nx][ny][nz][ndof]a =
  tabulate_3d nx ny nz (\i j k ->
    tabulate ndof (\dof ->
      if isOnBoundary (nx,ny,nz) {x=i,y=j,z=k} dof then
        zero
      else
        v[i,j,k,dof]))

def setBCtoInput [nx][ny][nz][ndof] 'a (input :[nx][ny][nz][ndof]a) (v :[nx][ny][nz][ndof]a) :[nx][ny][nz][ndof]a =
  tabulate_3d nx ny nz (\i j k ->
    tabulate ndof (\dof ->
      if isOnBoundary (nx,ny,nz) {x=i,y=j,z=k} dof then
        input[i,j,k,dof]
      else
        v[i,j,k,dof]))
