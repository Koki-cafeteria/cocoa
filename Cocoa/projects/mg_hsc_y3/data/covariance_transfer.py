# Python script to read covariance.dat, apply the specified transformations, and write to covariance2.dat

import numpy as np

# Define file names
input_file = "/home/tanida/GL/covariance/covariance.dat"
output_file = "covariance_transformed.dat"

# Read the file while skipping the first two comment lines
with open(input_file, 'r') as f:
    comments = [next(f) for _ in range(2)]  # First two lines are comments
    data = np.loadtxt(f)  # Load the remaining data as a NumPy array

# Verify data shape
if data.shape != (240, 240):
    raise ValueError(f"Expected data shape (240, 240), but got {data.shape}")

# Create a new empty array for transformed data
transformed_data = np.zeros_like(data)

# Apply the transformations as specified
transformed_data[60:150, 60:150] = data[0:90, 0:90]
transformed_data[0:60, 60:150] = data[90:150, 0:90]
transformed_data[150:240, 60:150] = data[150:240, 0:90]

transformed_data[60:150, 0:60] = data[0:90, 90:150]
transformed_data[0:60, 0:60] = data[90:150, 90:150]
transformed_data[150:240, 0:60] = data[150:240, 90:150]

transformed_data[60:150, 150:240] = data[0:90, 150:240]
transformed_data[0:60, 150:240] = data[90:150, 150:240]
transformed_data[150:240, 150:240] = data[150:240, 150:240]

# Write the transformed data to the output file, including the original comments
with open(output_file, 'w') as f:
    f.writelines(comments)  # Write the comment lines
    np.savetxt(f, transformed_data, fmt='%.6e')  # Save the transformed data

output_file

