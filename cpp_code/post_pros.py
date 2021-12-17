import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib

def main():
	"""! The post processing file that reads the desinty field exported by the topology optimization to produce a pretty picture.


	This program takes as input from the command line the number of elements along x and along y and uses the density_field.csv present
 	the current working directory to produce the density_field.png as a grayscale color map plot. This helps visualise the topology optimized
 	domain
	"""
	# Check if we have the required number of arguments
	if(len(sys.argv) < 2):
		print(f"Please Enter atleast 2 arguments to run this program")
		sys.exit(0)

	# Get the data from the csv file into a numpy array
	# Set plot attributes
	data = np.genfromtxt("density_field.csv",delimiter = ',')
	font = {'family' : 'normal',
	        'weight' : 'bold',
	        'size'   : 18}
	matplotlib.rc('font', **font)


	# Generate figure
	fig= plt.plot(figsize = (10,8))
	# Define x and y
	x = float(sys.argv[1])
	y = float(sys.argv[2])
	# Plot the color map in grayscale
	plt.imshow(data, cmap='gray',extent=[0,x,0,y], aspect=1)
	plt.title("Density Field")
	plt.savefig("density_field.png")
	plt.show()



if __name__ == '__main__':
	main()