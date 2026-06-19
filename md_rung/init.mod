variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5
units metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2
variable a equal 3.615
boundary p p p
lattice fcc ${a}
region box prism 0 2.0 0 3.0 0 4.0 0.0 0.0 0.0
create_box 1 box
create_atoms 1 box
mass 1 63.546
