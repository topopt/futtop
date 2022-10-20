-- n-dimensional map
def map_2d 'a 'x [n][m] (f: a -> x) (as: [n][m]a): [n][m]x = map (map f) as
def map_3d 'a 'x [n][m][l] (f: a -> x) (as: [n][m][l]a): [n][m][l]x = map (map_2d f) as
def map_4d 'a 'x [n][m][l][k] (f: a -> x) (as: [n][m][l][k]a): [n][m][l][k]x = map (map_3d f) as

-- n-dimensional map2
def map2_2d 'a 'b 'x [n][m] (f: a -> b -> x) (as: [n][m]a) (bs: [n][m]b): [n][m]x =
  map2 (map2 f) as bs
def map2_3d 'a 'b 'x [n][m][l] (f: a -> b -> x) (as: [n][m][l]a) (bs: [n][m][l]b): [n][m][l]x =
  map2 (map2_2d f) as bs
def map2_4d 'a 'b 'x [n][m][l][k] (f: a -> b -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b): [n][m][l][k]x =
  map2 (map2_3d f) as bs

-- n-dimensional map3
def map3_3d 'a 'b 'c 'x [n][m][l] (f: a -> b -> c -> x) (as: [n][m][l]a) (bs: [n][m][l]b) (cs: [n][m][l]c) : [n][m][l]x =
  map3 (map3 (map3 f)) as bs cs
def map3_4d 'a 'b 'c 'x [n][m][l][k] (f: a -> b -> c -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b) (cs: [n][m][l][k]c) : [n][m][l][k]x =
  map3 (map3_3d f) as bs cs

-- 4-dimensional map4
def map4_2d 'a 'b 'c 'd 'x [n][m] (f: a -> b -> c -> d -> x) (as: [n][m]a) (bs: [n][m]b) (cs: [n][m]c) (ds: [n][m]d): [n][m]x =
  map4 (map4 f) as bs cs ds
def map4_3d 'a 'b 'c 'd 'x [n][m][l] (f: a -> b -> c -> d -> x) (as: [n][m][l]a) (bs: [n][m][l]b) (cs: [n][m][l]c) (ds: [n][m][l]d): [n][m][l]x =
  map4 (map4_2d f) as bs cs ds
def map4_4d 'a 'b 'c 'd 'x [n][m][l][k] (f: a -> b -> c -> d -> x) (as: [n][m][l][k]a) (bs: [n][m][l][k]b) (cs: [n][m][l][k]c) (ds: [n][m][l][k]d): [n][m][l][k]x =
  map4 (map4_3d f) as bs cs ds


def replicate_2d a b v = replicate a (replicate b v)
def replicate_3d a b c v = replicate a (replicate_2d b c v)
def replicate_4d a b c d v = replicate a (replicate_3d b c d v)

def vecmul [n][m] 'a (add: a -> a -> a) (mul: a -> a -> a) (zero: a) (A: [n][m]a) (x: [m]a) : [n]a =
  map (\A_row -> reduce add zero (map2 mul A_row x)) A

def vecmul_f64 = vecmul (+) (*) 0f64
