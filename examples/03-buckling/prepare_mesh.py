import igl

fname = "column2"
v, f = igl.read_triangle_mesh(fname+".obj")

file = open(fname+".gts", "w")

# print statistics
file.write(str(len(v))+" "+str(len(f))+" \n")

for ve in v:
	file.write(str(ve[0])+" "+str(ve[1])+" "+str(ve[2])+"\n")

for fe in f:
	file.write(str(fe[0])+" "+str(fe[1])+" "+str(fe[2])+"\n")

file.close()