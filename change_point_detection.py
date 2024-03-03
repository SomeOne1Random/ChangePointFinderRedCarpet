import matplotlib.pyplot as plt
import pandas as pd
import ruptures as rpt
from multiprocessing import Pool
import numpy as np
import os
from scipy import stats

def detect_change_points(data_chunk):
    algo = rpt.KernelCPD(kernel="linear").fit(data_chunk)
    change_points_with_penalty = algo.predict(pen=700)
    if len(change_points_with_penalty) > 1:
        return algo.predict(n_bkps=2)
    return change_points_with_penalty

def process_segment(data_chunk):
    result = detect_change_points(data_chunk)
    return result

def calculate_region_features(data, start, end):
    region_data = data[:, start:end]
    mean = np.mean(region_data)
    variance = np.var(region_data)
    return mean, variance

def calculate_region_t_value(data, region1, region2):
    start1, end1 = region1
    start2, end2 = region2

    data_region1 = data[:, start1:end1].flatten()
    data_region2 = data[:, start2:end2].flatten()

    # Perform independent t-test
    t_stat, p_value = stats.ttest_ind(data_region1, data_region2, equal_var=False)
    return t_stat, p_value

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
def write_regions_to_file(regions, file_path):
    with open(file_path, 'w') as file:
        for region_name, (start, end) in regions.items():
            file.write(f"{region_name}: {start} - {end}\n")

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

    all_change_points = [0]
    for segment_result in results:
        for bkpt in segment_result[:-1]:
            all_change_points.append(bkpt)
    all_change_points.append(len(redcarpet_npy[0]))
    all_change_points = sorted(list(set(all_change_points)))

    regions = {}
    for i in range(len(all_change_points) - 1):
        region_name = f"Region {i+1}"
        regions[region_name] = (all_change_points[i], all_change_points[i+1])

    plt.figure(figsize=(20, 10))
    for region_name, (start, end) in regions.items():
        plt.axvspan(start, end, color=np.random.rand(3,), alpha=0.3, label=region_name)
    plt.plot(redcarpet_npy[0, :], lw=1)
    plt.title('Data Visualization with Regions')
    plt.legend()
    plt.show()

    directory = os.path.dirname(file_path)
    output_file_path_change_points = os.path.join(directory, 'ChangePoints.txt')
    output_file_path_regions = os.path.join(directory, 'regions.txt')
    write_regions_to_file(regions, output_file_path_regions)

    # New steps for comparing regions
    write_change_points_to_file(results, output_file_path_change_points)

    # Calculate and write region T-values to file
    with open(output_file_path_regions, 'a') as file:
        file.write("\nRegion T-values:\n")
        for region_name1, (start1, end1) in regions.items():
            for region_name2, (start2, end2) in regions.items():
                if region_name1 != region_name2:
                    t_stat, p_value = calculate_region_t_value(redcarpet_npy, (start1, end1), (start2, end2))
                    file.write(f"{region_name1} vs {region_name2}: T-value = {t_stat}, P-value = {p_value}\n")


if __name__ == "__main__":
    main()
