# Importing necessary libraries
import matplotlib.pyplot as plt  # for plotting
import pandas as pd  # for data manipulation
import ruptures as rpt  # for change point detection
from multiprocessing import Pool  # for parallel processing
import numpy as np  # for numerical operations
import os  # for file and directory operations
'''
Finds Changepoints in a RedCarpet genome
Paramters are pen and n_bkps 
Additionally, the amount of cores used can be adjusted with num_chunks
'''

# Function to detect change points in a given chunk of data
def detect_change_points(data_chunk):
    # Initializing the change point detection algorithm with a linear kernel
    algo = rpt.KernelCPD(kernel="linear").fit(data_chunk)

    # Detecting change points with a specified penalty value
    change_points_with_penalty = algo.predict(pen=1000)  # pen=1000 is an adjustable parameter

    # Using a different method for change point detection if the penalty method finds more than one change point
    if len(change_points_with_penalty) > 1:
        return algo.predict(n_bkps=3)  # n_bkps=3 is an adjustable parameter

    # Returning the results from the penalty method if no additional change points are found
    return change_points_with_penalty

# Function to process each data segment for change point detection
def process_segment(data_chunk):
    # Detecting change points in the provided data chunk
    result = detect_change_points(data_chunk)
    return result

# Function to write the detected change points to a file
def write_change_points_to_file(results, file_path):
    # Aggregating all change points from different segments, excluding the last point in each segment
    all_change_points = set()
    for segment_result in results:
        all_change_points.update(segment_result[:-1])

    # Sorting the change points for readability
    sorted_change_points = sorted(all_change_points)

    # Writing the sorted change points to the specified file
    with open(file_path, 'w') as file:
        file.write("Detected Change Points: " + ", ".join(map(str, sorted_change_points)) + "\n")

# Main function to execute the script
def main():
    # Path to the input data file
    file_path = '/Users/srujanyamali/Downloads/Microbal Stuff/Microbal Datasets/TB_bk_all_Redcarpet_report.txt' #file path

    # Reading the data from the file using pandas
    redcarpet = pd.read_csv(file_path, sep='\t')
    # Converting the dataframe to a numpy array and ensuring the data type is float
    redcarpet_npy = redcarpet.to_numpy().astype('float')

    # Number of chunks to split the data into
    num_chunks = 10
    # Splitting the data array into specified number of chunks
    data_chunks = np.array_split(redcarpet_npy, num_chunks, axis=1)

    # Creating a pool of processes for parallel execution
    with Pool(processes=num_chunks) as pool:
        # Processing each data chunk in parallel and storing the results
        results = pool.map(process_segment, data_chunks)

    # Setting up the plot for visualization
    plt.figure(figsize=(20, 20))
    # Creating a heatmap of the data
    plt.imshow(redcarpet_npy, cmap='hot', interpolation='nearest')

    # Adding lines to the plot for each detected change point
    for segment_result in results:
        for bkpt in segment_result[:-1]:
            plt.axvline(x=bkpt, color='cyan', linestyle='-')
            plt.axhline(y=bkpt, color='cyan', linestyle='-')

    # Setting labels and title for the plot
    plt.title('Visualization of Protein Changepoints with Detected Change Points')
    plt.ylabel('Protein [ordered]')
    plt.xlabel('Protein [ordered]')
    # Displaying the plot
    plt.show()

    # Getting the directory of the input file
    directory = os.path.dirname(file_path)
    # Constructing the path for the output file
    output_file_path = os.path.join(directory, 'change_points.txt')

    # Writing the detected change points to the output file
    write_change_points_to_file(results, output_file_path)

# Condition to check if the script is being run as the main program
if __name__ == "__main__":
    main()  # Executing the main function
