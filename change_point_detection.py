# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 22:28:37 2021
@author: Arnav Aryaraj
#made with <3 in Atlanta, GA and Phila, PA
"""
#Carpet Cleaner: Finds the useful "stains" (squares) in the carpet
#%%imports
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import find_peaks
start_time = time.time()
#%% PARAMETERS; editable terms to be used in the rest of the script.

#IMPORTS
redcarpet = pd.read_csv('/Users/srujanyamali/Downloads/Microbal Stuff/Microbal Datasets/recombination_report_TW20_intersection_all.txt',sep='\t') #import redcarpet

#PARAMETERS [these values are somewhat arbitrary, and can be changed. Larger windows might be better at picking up larger recombinations, lower height_peaks might pick up older recombinations, etc.]
window = 100 #how many proteins on either side to account for when performing changepoint analysis
height_peaks = 100000 #look for peaks above this threshold of log likelihood. If you change the window, you will probably have to change this too.
dist = 100000 #if peaks are too close to each other (i.e., the distance between two local peaks < 'dist'), keep only the higher peak. ...
          #... This is designed to remove multiple peaks where changepoints likelihood is very high

#%% PREPARE CARPET
#import carpet file, and extend the array by add a little bit of the end to the beginning and vice versa.
#bacterial genomes are circular, so this will allow the first protein to be analyzed by creating a window on the LHS comprised of the last 20 proteins in the file.

#Parameters: 'redcarpet', 'window'

start_time = time.time()
redcarpet_npy = redcarpet.to_numpy().astype('float') #convert redcarpet to numpy array
rc_length = len(redcarpet_npy)
redcarpet_npyext = np.zeros((rc_length+window+window, rc_length)) #create an array with a little bit of buffer on top and bottom.
redcarpet_npyext[0:window,:] = redcarpet_npy[rc_length-window:rc_length,:] #fill the top of the extended array
redcarpet_npyext[len(redcarpet_npyext)-window:len(redcarpet_npyext),:] = redcarpet_npy[0:window,:] #fill the bottom of the extended array
redcarpet_npyext[window:window+rc_length,:] = redcarpet_npy #fill the middle of the extended array
#plt.imshow(redcarpet_npyext, cmap='hot', interpolation='nearest')  #if you want to see what the extended array looks like
print("Done Preparing Carpet for Analysis:--- %s seconds ---" % (time.time() - start_time))
#%%   EXTENDED VISUALIZATION -- UNNECESSARY
# as the visual would indicate the specific choice of j is essentially irrelevant
# the changepoint values are identical for any specific value of j chosen.
#so 'fast_score' only calculates all values of i when j = 0.
# this section just creates visuals that prove this point.
"""
score = np.zeros((rc_length,rc_length))
for i in tqdm(np.arange(rc_length)):
    for j in np.arange(rc_length):
        temp = redcarpet_npyext[:,i]
        k = i + window
        score_early = redcarpet_npyext[k-window:k]
        score_late = redcarpet_npyext[k:k+window]
        score_tot = redcarpet_npyext[k-window:k+window]
        score[j,i] = ((np.sum(score_late))*np.log((np.average(score_late))))+((np.sum(score_early))*np.log((np.average(score_early))))-((np.sum(score_tot))*np.log((np.average(score_tot))))
plt.figure(0,figsize=(20, 20))
plt.imshow(score, cmap='hot', interpolation='nearest')
plt.xlabel('Protein Assesed for Changepoint')
plt.ylabel('Column used to assess Protein Chanepoint')
plt.title('Changepoint Values: Redundancy in Columns for Changepoint analysis')
plt.colorbar()
"""
#%% SCORE
#assign likelihoods of changepoint by looking at the +/- window space for a changepoint.

#Parameters: 'window'

fast_score = np.zeros((1,rc_length)) #calculate changepoint score with one column (irrelevant which one, so first one; j=0)
for i in tqdm(np.arange(rc_length)):
    for j in np.arange(1): #chose first column: j=0
        temp = redcarpet_npyext[:,i] #consider the ith slice of the redcarpet
        k = i + window #start a window away, so that you can calculate score_early (i.e., the -20th term; this is why the buffer was added).
        score_early = redcarpet_npyext[k-window:k] #the prior window of points
        score_late = redcarpet_npyext[k:k+window] #the next window of points
        score_tot = redcarpet_npyext[k-window:k+window] #both the prior and next windows
        fast_score[j,i] = ((np.sum(score_late))*np.log((np.average(score_late))))+((np.sum(score_early))*np.log((np.average(score_early))))-((np.sum(score_tot))*np.log((np.average(score_tot))))
        # ^ changepoint analysis; https://repository.upenn.edu/cgi/viewcontent.cgi?article=1520&amp;context=physics_papers
print("Changepoint Likelihood Scores Calculated:--- %s seconds ---" % (time.time() - start_time))
#%%
plt.figure(1,figsize=(20, 6))
plt.plot(fast_score[0,:])
plt.xlabel('Protein [ordered]')
plt.ylabel('Changepoint Log Likelihood [au]')
plt.title('Likelihood of Changepoint within Genome')

#extend the score (take a bit of the LHS and add to the RHS and vice versa. This is important for finding peaks, so ...
# ...that peaks near the beginning/end are accounted for.
ext_fast_score = np.zeros((1,rc_length+window+window))
ext_fast_score[0,0:window] = fast_score[0,rc_length-window:rc_length]
ext_fast_score[0,len(redcarpet_npyext)-window:len(redcarpet_npyext)] = fast_score[0,0:window]
ext_fast_score[0,window:window+rc_length] = fast_score[0,:]
#%% FIND CHANGEPOINT PEAKS

# Parameters: 'height_peaks'

# Smooth the likelihood scores to reduce noise
smoothed_score = np.convolve(ext_fast_score[0, :], np.ones(5) / 5, mode='same')

# Find peaks in the smoothed scores
peaks, _ = find_peaks(smoothed_score, height=height_peaks)

# Refine the peaks using peak width and prominence
peaks, _ = find_peaks(smoothed_score, height=height_peaks, width=10, prominence=10)

# Move peaks back to the true positions
peaks = peaks - window

# Remove peaks that are too close to the beginning or end
peaks = [p for p in peaks if p >= 0 and p < rc_length]

print("Changepoint Likelihood Scores Calculated:--- %s seconds ---" % (time.time() - start_time))

#%%POLISH CHANGEPOINT
# This section exists becausse sometimes, the peaks are not perfect, and multiple local maxima exist within a small space
# i.e., peak 2942, 2945 and 2947are all local maxima. Usually, this is probably because one of the three is the true most likely changepoint.
# noise caused the other two points to be maxima. So, this section will remove those other two points.

# Parameters: 'dist'

peaks_rm = [] #peaks to remove
numerical_list = list(peaks)
numerical_list.sort() #sort peaks in order of appearance
for i in np.arange(len(numerical_list)-1):
    if abs(numerical_list[i] - numerical_list[i+1]) < dist:  #if two peaks are close
        if fast_score[0,:][numerical_list[i]] > fast_score[0,:][numerical_list[i+1]]: #see which peak is higher
            peaks_rm.append(numerical_list[i+1]) #remove the lower peak
        else:
            peaks_rm.append(numerical_list[i]) #remove the lower peak
peaks_rm = list(set(peaks_rm))
peaks_v2 = list(set(peaks).difference(peaks_rm)) #actually remove the lower local peaks from the list of peaks.
peaks_v2.sort()
print("Changepoint Likelihood Scores Polished:--- %s seconds ---" % (time.time() - start_time))

#%%VISUALIZATION

plt.figure(2, figsize=(20, 20))
plt.imshow(redcarpet, cmap='hot', interpolation='nearest')
for i in np.arange(0,len(peaks_v2)):
    plt.plot((peaks_v2[i])*np.ones(len(redcarpet)),np.arange(0,len(redcarpet)),'b',lw=2)
    plt.plot(np.arange(0,len(redcarpet)),(peaks_v2[i])*np.ones(len(redcarpet)),'b',lw=2)
plt.title('Visualization of Protein Changepoints')
plt.ylabel('Protein [ordered]')
plt.xlabel('Protein [ordered]')


plt.figure(3, figsize=(20, 6))
plt.plot(fast_score[0,:], label='Log Likelihood')
plt.plot(np.arange(0,len(fast_score[0,:])),height_peaks*np.ones(len(fast_score[0,:])), label = 'Threshold for Peaks')
plt.plot(peaks, fast_score[0,:][peaks],'o', label='Removed Local peaks',color='maroon')
plt.plot(peaks_v2, fast_score[0,:][peaks_v2],'o', label = 'Kept local peaks',color='palegreen')
plt.title('Finding Peaks Within theLikelihood of Changepoint within Genome')
plt.xlabel('Protein [ordered]')
plt.ylabel('Changepoint Log Likelihood [au]')
plt.legend()
plt.show()
print("Finished analysis:--- %s seconds ---" % (time.time() - start_time))
print ('')
print ('')
print ('All peaks are proteins: ' + str(list(peaks_v2)))
