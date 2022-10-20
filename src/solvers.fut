import "applyStiffnessMatrix"
import "multigrid"
import "utility"

-- added kahan summation from the futhark examples
type kahan = (f32, f32)
def kahan_add ((s1, c1) : kahan) ((s2, c2) : kahan) : kahan =
  let s = s1 + s2
  let d  = if f32.abs s1 >= f32.abs s2 then (s1 - s) + s2
           else (s2 - s) + s1
  let c = (c1 + c2) + d
  in (s, c)
def kahan_sum (xs: []f32) : f32 =
  let (s,c) = reduce kahan_add (0,0) (map (\x -> (x,0)) xs)
  in s + c

def innerProduct [n][m][l][k] (a :[n][m][l][k]f32) (b :[n][m][l][k]f32) :f32 =
  map2_4d (*) a b
  |> flatten_4d
  |> kahan_sum

def norm [n][m][l][k] (a :[n][m][l][k]f32) :f32 =
  map_4d (\x -> x*x) a
  |> flatten_4d
  |> kahan_sum
  |> f32.sqrt

def cgSolveMG [nelx][nely][nelz][nx][ny][nz] (x :[nelx][nely][nelz]f32) (mgData :multigridData) (b :[nx][ny][nz][3]f32) (u :[nx][ny][nz][3]f32) =
  let tol     = 1e-5
  let maxIt   = 200
  let bnorm   = norm b
  let zero_4d = replicate_4d nx ny nz 3 0f32

  let r = applyStiffnessMatrixSingle x u
    |> map2_4d (\bb rr -> bb - rr) b

  let  (u, _, _, relres, it, _) =
    loop (uold, rold, pold, res, its, rhoold) = (u, r, zero_4d, 1f32, 0i32, 1f32) while (res > tol && its < maxIt) do
      let z      = vcycle_l0 mgData x rold
      let rho    = innerProduct rold z
      let beta   = rho / rhoold
      let p      = map2_4d (\pp zz -> beta * pp + zz) pold z
      let q      = applyStiffnessMatrixSingle x p
      let alpha  = rho / (innerProduct p q)
      let u      = map2_4d (\uu pp -> uu + alpha * pp) uold p
      let r      = map2_4d (\rr qq -> rr - alpha * qq) rold q
      let relres = (norm r) / bnorm
    in (u, r, p, relres, its+1, rho)
  in (u, relres, it)
