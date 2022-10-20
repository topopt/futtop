import "indexUtilities"
import "keConstants"
import "material"

def getElementUnscaledCompliance [nx][ny][nz] (elementIndex :index) (u :[nx][ny][nz][3]f32) :f32 =
  let ulocal = readLocalState u elementIndex |> map f64.f32
  let Ku = keprod 1 (copy ulocal)
  in f32.f64 (reduce (+) 0 (map2 (*) Ku ulocal))

entry getComplianceSensitivity [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (u :[nx][ny][nz][3]f32)  :(f32, [nelx][nely][nelz]f32) =
  let unscaledCompliance = tabulate_3d nelx nely nelz (\i j k -> getElementUnscaledCompliance {x=i,y=j,z=k} u)
    |> flatten
    |> flatten_to (nelx*nely*nelz)
  let E    = f32.f64 E
  let Emin = f32.f64 Emin
  let xflat = x
    |> flatten
    |> flatten_to (nelx*nely*nelz)
  let c     = map2 (\rho uku -> (Emin + (E-Emin)*rho*rho*rho)*uku) xflat unscaledCompliance |> reduce (+) 0
  let csens = map2 (\rho uku -> (      (-3)*(E-Emin)*rho*rho)*uku) xflat unscaledCompliance |> unflatten_3d nelx nely nelz
  in (c, csens)

-- ==
-- entry: getComplianceSensitivity
-- nobench input @../testData/getComplianceSens.txt output @../testData/getComplianceSensOut.txt
-- compiled random input { [64][64][64]f32 [65][65][65][3]f32 }
-- compiled random input { [128][128][128]f32 [129][129][129][3]f32 }
