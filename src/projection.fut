def projectToFiner [nx][ny][nz] (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  let nxf = ((nx-1)*2)+1
  let nyf = ((ny-1)*2)+1
  let nzf = ((nz-1)*2)+1
  in tabulate_3d nxf nyf nzf (\i j k ->
      let ic1 =  (i)     / 2
      let ic2 =  (i + 1) / 2
      let jc1 =  (j)     / 2
      let jc2 =  (j + 1) / 2
      let kc1 =  (k)     / 2
      let kc2 =  (k + 1) / 2
      in #[unsafe] ([u[ic1,jc1,kc1],u[ic1,jc1,kc2],
                     u[ic1,jc2,kc1],u[ic1,jc2,kc2],
                     u[ic2,jc1,kc1],u[ic2,jc1,kc2],
                     u[ic2,jc2,kc1],u[ic2,jc2,kc2]])
          |> transpose
          |> map (f64.sum)
          |> map (*0.125))

def loadValue 't [nx][ny][nz] (u :[nx][ny][nz]t) (zero :t) (idx :(i64,i64,i64)) =
  if (idx.0 < 0 || idx.1 < 0 || idx.2 < 0 || idx.0 >= nx || idx.1 >= ny || idx.2 >= nz)
  then zero
  else #[unsafe] (u[idx.0,idx.1,idx.2])

def getWeight (idx :(i64,i64,i64)) :f64 =
  let dist = i64.abs(idx.0) + i64.abs(idx.1) + i64.abs(idx.2)
  in match dist
    case 0 -> 1.0
    case 1 -> 0.5
    case 2 -> 0.25
    case 3 -> 0.125
    case _ -> 0.0

def projectToCoarser [nx][ny][nz] (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  let nxc = ((nx-1)/2)+1
  let nyc = ((ny-1)/2)+1
  let nzc = ((nz-1)/2)+1
  in tabulate_3d nxc nyc nzc (\i j k ->
      let icenter = i * 2
      let jcenter = j * 2
      let kcenter = k * 2
      let offsets =
      [(-1,-1,-1), (-1,-1, 0), (-1,-1, 1),
       (-1, 0,-1), (-1, 0, 0), (-1, 0, 1),
       (-1, 1,-1), (-1, 1, 0), (-1, 1, 1),
       ( 0,-1,-1), ( 0,-1, 0), ( 0,-1, 1),
       ( 0, 0,-1), ( 0, 0, 0), ( 0, 0, 1),
       ( 0, 1,-1), ( 0, 1, 0), ( 0, 1, 1),
       ( 1,-1,-1), ( 1,-1, 0), ( 1,-1, 1),
       ( 1, 0,-1), ( 1, 0, 0), ( 1, 0, 1),
       ( 1, 1,-1), ( 1, 1, 0), ( 1, 1, 1)]
      let weights = offsets |> map getWeight
      in offsets
        |> map(\p -> (p.0+icenter,p.1+jcenter,p.2+kcenter))
        |> map(loadValue u [0,0,0])
        |> map2(\w x -> map (\y -> w*y) x) weights
        |> transpose
        |> map (f64.sum))


def projectToFinerSingle [nx][ny][nz] (u :[nx][ny][nz][3]f32) :[][][][3]f32 =
  let nxf = ((nx-1)*2)+1
  let nyf = ((ny-1)*2)+1
  let nzf = ((nz-1)*2)+1
  in tabulate_3d nxf nyf nzf (\i j k ->
      let ic1 =  (i)     / 2
      let ic2 =  (i + 1) / 2
      let jc1 =  (j)     / 2
      let jc2 =  (j + 1) / 2
      let kc1 =  (k)     / 2
      let kc2 =  (k + 1) / 2
      in #[unsafe] ([u[ic1,jc1,kc1],u[ic1,jc1,kc2],
                     u[ic1,jc2,kc1],u[ic1,jc2,kc2],
                     u[ic2,jc1,kc1],u[ic2,jc1,kc2],
                     u[ic2,jc2,kc1],u[ic2,jc2,kc2]])
          |> transpose
          |> map (f32.sum)
          |> map (*0.125))


def projectToCoarserSingle [nx][ny][nz] (u :[nx][ny][nz][3]f32) :[][][][3]f32 =
  let nxc = ((nx-1)/2)+1
  let nyc = ((ny-1)/2)+1
  let nzc = ((nz-1)/2)+1
  in tabulate_3d nxc nyc nzc (\i j k ->
      let icenter = i * 2
      let jcenter = j * 2
      let kcenter = k * 2
      let offsets =
      [(-1,-1,-1), (-1,-1, 0), (-1,-1, 1),
       (-1, 0,-1), (-1, 0, 0), (-1, 0, 1),
       (-1, 1,-1), (-1, 1, 0), (-1, 1, 1),
       ( 0,-1,-1), ( 0,-1, 0), ( 0,-1, 1),
       ( 0, 0,-1), ( 0, 0, 0), ( 0, 0, 1),
       ( 0, 1,-1), ( 0, 1, 0), ( 0, 1, 1),
       ( 1,-1,-1), ( 1,-1, 0), ( 1,-1, 1),
       ( 1, 0,-1), ( 1, 0, 0), ( 1, 0, 1),
       ( 1, 1,-1), ( 1, 1, 0), ( 1, 1, 1)]
      let weights = offsets |> map getWeight |> map f32.f64
      in offsets
        |> map(\p -> (p.0+icenter,p.1+jcenter,p.2+kcenter))
        |> map(loadValue u [0,0,0])
        |> map2(\w x -> map (\y -> w*y) x) weights
        |> transpose
        |> map (f32.sum))


entry projectToFiner_test [nx][ny][nz] (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  projectToFiner u

entry projectToCoarser_test [nx][ny][nz] (u :[nx][ny][nz][3]f64) :[][][][3]f64 =
  projectToCoarser u

-- ==
-- entry: projectToFiner_test
-- compiled input @../testData/projectToCoarser1Out.txt output @../testData/projectToFiner1Out.txt
-- compiled input @../testData/projectToCoarser2Out.txt output @../testData/projectToFiner2Out.txt
-- compiled input @../testData/projectToCoarser3Out.txt output @../testData/projectToFiner3Out.txt
-- compiled random input { [33][33][33][3]f64 } auto output
-- compiled random input { [65][65][65][3]f64 }
-- compiled random input { [129][129][129][3]f64 }

-- ==
-- entry: projectToCoarser_test
-- compiled input @../testData/projectToCoarser1.txt output @../testData/projectToCoarser1Out.txt
-- compiled input @../testData/projectToCoarser2.txt output @../testData/projectToCoarser2Out.txt
-- compiled input @../testData/projectToCoarser3.txt output @../testData/projectToCoarser3Out.txt
-- compiled random input { [65][65][65][3]f64 } auto output
-- compiled random input { [129][129][129][3]f64 }
-- compiled random input { [257][257][257][3]f64 }
