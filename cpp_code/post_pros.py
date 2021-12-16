import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib

def main():
	if(len(sys.argv) < 2):
		print(f"Please Enter atleast 2 arguments to run this program")
		sys.exit(0)

	data = np.genfromtxt("density_field.csv",delimiter = ',')
	font = {'family' : 'normal',
	        'weight' : 'bold',
	        'size'   : 18}
	matplotlib.rc('font', **font)



	fig= plt.plot(figsize = (10,8))
	print(f"{sys.argv[1]}, {sys.argv[2]}")
	x = float(sys.argv[1])
	y = float(sys.argv[2])
	plt.imshow(data, cmap='gray',extent=[0,x,0,y], aspect=1)

	plt.title("Density Field")
	plt.savefig("density_field.png")
	plt.show()



if __name__ == '__main__':
	main()