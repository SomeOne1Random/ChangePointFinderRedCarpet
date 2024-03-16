import matplotlib.pyplot as plt
import pandas as pd
import ruptures as rpt
from multiprocessing import Pool
import numpy as np
import os

def detect_change_points(data_chunk):
    algo = rpt.KernelCPD(kernel="linear").fit(data_chunk)

    # Check for change points with a penalty value of 700
    change_points_with_penalty = algo.predict(pen=700)

    # If change points are found with the penalty method, use n_bkps for actual detection
    if len(change_points_with_penalty) > 1:
        return algo.predict(n_bkps=2)

        # If no change points are found with the penalty method, return the penalty results
    return change_points_with_penalty


def process_segment(data_chunk):
    result = detect_change_points(data_chunk)
    return result

def write_change_points_to_file(results, file_path):
    # Flatten the results into a single list
    all_change_points = set()
    for segment_result in results:
        all_change_points.update(segment_result[:-1])

        # Sort the change points for better readability
    sorted_change_points = sorted(all_change_points)

    # Write the change points to the file
    with open(file_path, 'w') as file:
        file.write("Detected Change Points: " + ", ".join(map(str, sorted_change_points)) + "\n")

def main():
    file_path = '/Users/srujanyamali/Downloads/Microbial Stuff/Microbial Datasets/GCA_000027045.1_ASM2704v1_modified_Redcarpet_report.txt'
    redcarpet = pd.read_csv(file_path, sep='\t')
    redcarpet_npy = redcarpet.to_numpy().astype('float')

    num_chunks = 10
    data_chunks = np.array_split(redcarpet_npy, num_chunks, axis=1)

    with Pool(processes=num_chunks) as pool:
        results = pool.map(process_segment, data_chunks)

    plt.figure(figsize=(20, 20))
    plt.imshow(redcarpet_npy, cmap='hot', interpolation='nearest')

    for segment_result in results:
        for bkpt in segment_result[:-1]:
            plt.axvline(x=bkpt, color='cyan', linestyle='-')
            plt.axhline(y=bkpt, color='cyan', linestyle='-')

    plt.title('Visualization of Protein Changepoints with Detected Change Points')
    plt.ylabel('Protein [ordered]')
    plt.xlabel('Protein [ordered]')
    plt.show()

    # Extract the first row (or any specific row or column) for plotting
    plot_data = redcarpet_npy[0, :]  # This gets the first row. Adjust the indices as needed.

    plt.figure(figsize=(40, 15))
    plt.plot(plot_data, lw=1)  # Plot the line

    for segment_result in results:
        for bkpt in segment_result[:-1]:
            if bkpt < len(plot_data):
                plt.axvline(x=bkpt, color='orange', linestyle='-',
                            linewidth=2)  # Draw a vertical line at each change point

    plt.title('Visualization of Change Points')
    plt.ylabel('Value')
    plt.xlabel('Index')
    plt.show()

    directory = os.path.dirname(file_path)

    # Create the output file path in the same directory
    output_file_path = os.path.join(directory, 'change_points.txt')

    # Write the change points to a file
    write_change_points_to_file(results, output_file_path)

if __name__ == "__main__":
    main()