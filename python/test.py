from dolfin import *
from use_meshfunction import get_dofs0
from use_indexsets import get_dofs1

class MovingCircle(Expression):
  def __init__(self, cx, cy, r, dt):
    Expression.__init__(self)
    self.cx = cx
    self.cy = cy
    self.r = r
    self.dt = dt

  def update(self, motion_type):
    if motion_type == "rotation":  
      self.cx, self.cy = self.cx - self.cy*self.dt, self.cy + self.cx*self.dt
    
    if motion_type == "translation":
      self.cx += self.dt

  def eval(self, values, x):
    values[0] = self.r - sqrt((x[0] - self.cx)**2 + (x[1] - self.cy)**2) 

#------------------------------------------------------------------------------

def compare():
  '''Compare get_dofs0 and get_dofs1. Should return empty set.'''
  mesh = RectangleMesh(-2, -2, 2, 2, 100, 100)
  V = FunctionSpace(mesh, "CG", 1)
  size = V.dim()

  circle = MovingCircle(-1, 0, 0.5, 0.1)

  u = Function(V)
  circle.update("translate")
  u.interpolate(circle)

  dofs0 = get_dofs0(u, [-0.1, 0.1], 3)

  dofs1 = get_dofs1(u, [-0.1, 0.1], 3)


#------------------------------------------------------------------------------

def time_dofs(method0, method1, num_cells, width, motion_type):
  '''Time the get_dofs methods.'''
  mesh = RectangleMesh(-2, -2, 2, 2, num_cells, num_cells)
  mesh.init(2, 2)
  V = FunctionSpace(mesh, "CG", 1)
  size = V.dim()

  circle = MovingCircle(-1, 0, 0.5, 0.1)

  u = Function(V)

  for step in range(10):
    circle.update(motion_type)
    u.interpolate(circle)

    dofs0 = method0(u, [-0.1, 0.1], width)
    dofs1 = method1(u, [-0.1, 0.1], width)


  print "Dofs differ on set with measure", set(dofs1) - (set(dofs0)) # 0 desired

  print "width %d, n_cells %d cells, time %g vs %f in dofs of which" % \
  (2*width, mesh.num_cells(), timing("Getting dofs0"), timing("Getting dofs1"))
  
  print "\t\t first search %g vs %g" %\
  (timing("First search0"), timing("First search1"))

  print "\t\t band build %g vs %g" %\
  (timing("Band build0"), timing("Band build1"))

  print "\t\t connectivity %g vs %g" %\
  (timing("Connectivity0"), timing("Connectivity1"))

  print "\t\t extract dofs %g vs %g" %\
  (timing("Extract dofs0"), timing("Extract dofs1"))
#------------------------------------------------------------------------------

for motion_type in ["rotation", "translation"]:
  print motion_type
  for num_cells in [75, 125, 200]:
    for width in [0, 1, 2, 3]:
        time_dofs(method0=get_dofs0, method1=get_dofs1, num_cells=num_cells,\
                  width=width, motion_type=motion_type)
