import "indexUtilities"
import "keConstants"
import "utility"

type isoCoordinates = {xi: f64, eta: f64, zeta: f64}

def E    :f64 = 1
def Emin :f64 = 1e-6

def E_f32    :f32 = f32.f64 E
def Emin_f32 :f32 = f32.f64 Emin

def getYoungsModule (x :f64) = (Emin + (E-Emin)*x*x*x)
def getYoungsModuleDerivative (x :f64) = 3*(E-Emin)*x*x

def nu     :f64 = 0.3
def lambda :f64 = nu*E / ((1+nu)*(1-2*nu))
def mu     :f64 = E / (2*(1+nu))

def getElementYoungsModule [nelx][nely][nelz] (x :[nelx][nely][nelz]f32) (elementIndex :index) :f64 =
  if (indexIsInside (nelx,nely,nelz) elementIndex) then
    let xloc :f64 = getDensityUnsafe x elementIndex
    in getYoungsModule xloc
  else 0

def keprod (s :f64) (f :*[24]f64) :*[24]f64 =
  map (\keRow -> reduce (+) 0 (map2(\ff kr -> ff*kr*s) f keRow)) keconst
