import argparse
import scipy.special
import csv

parser = argparse.ArgumentParser(description="""Evaluate clustering results. This scripts will return the 
                                                precision, recall, F1-score and ARI of the provided clustering result""")

parser.add_argument("--clustered", 
                    required=True, 
                    type=str,
                    help="path to the .csv file with the clustering result")

parser.add_argument("--goldstandard", 
                    required=True,
                    type=str,
                    help="path to the .csv file with the gold standard")

args = vars(parser.parse_args())

# Get paths to clustering result and gold standard
clustered_file = args["clustered"]
gold_standard_file = args["goldstandard"]
# Save all the print in a list to save in a file later
results = []
results.append("{}{}\n".format("Clustering results file:",clustered_file))
results.append("{}{}\n".format("gold standard file:",gold_standard_file))

#print("\nStarting evaluate.py...")
#print("Clustering results file:", clustered_file)
#print("gold standard file:", gold_standard_file)

# Get the number of clusters from the gold standard
#---------------------------------------------------------
gold_standard_n_clusters = 0

all_gold_standard_clusters_list = []

with open(gold_standard_file) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        all_gold_standard_clusters_list.append(row[1])
        
gold_standard_clusters_list = list(set(all_gold_standard_clusters_list))
gold_standard_n_clusters = len(gold_standard_clusters_list)

results.append("{}{}\n".format("\nNumber of clusters available in the gold standard: ",gold_standard_n_clusters))

#print("\nNumber of clusters available in the gold standard: ", gold_standard_n_clusters)


# Get the gold standard
#----------------------------
gold_standard_clusters = [[] for x in range(gold_standard_n_clusters)]

gold_standard_count = 0

with open(gold_standard_file) as contig_clusters:
    readCSV = csv.reader(contig_clusters, delimiter=',')
    for row in readCSV:
        gold_standard_count += 1
        contig = row[0]
        bin_num = gold_standard_clusters_list.index(row[1])
        gold_standard_clusters[bin_num].append(contig)

results.append("{}{}\n".format("Number of objects available in the gold standard: ",gold_standard_count))

#print("Number of objects available in the gold standard: ", gold_standard_count)

# Get the number of clusters from the initial clustering result
#---------------------------------------------------------
n_clusters = 0

all_clusters_list = []

with open(clustered_file) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        all_clusters_list.append(row[1])
        
clusters_list = list(set(all_clusters_list))
n_clusters = len(clusters_list)

results.append("{}{}\n".format("Number of clusters available in the clustering result: ",n_clusters))
#print("Number of clusters available in the clustering result: ", n_clusters)


# Get initial clustering result
#----------------------------
clusters = [[] for x in range(n_clusters)]

clustered_count = 0
clustered_objects = []

with open(clustered_file) as contig_clusters:
    readCSV = csv.reader(contig_clusters, delimiter=',')
    for row in readCSV:
        clustered_count += 1
        contig = row[0]
        bin_num = clusters_list.index(row[1])
        clusters[bin_num].append(contig)
        clustered_objects.append(contig)

results.append("{}{}\n".format("Number of objects available in the clustering result: ",len(clustered_objects)))
#print("Number of objects available in the clustering result: ", len(clustered_objects))

# Functions to determine precision, recall, F1-score and ARI
#------------------------------------------------------------

# Get precicion
def getPrecision(mat, k, s, total):
    sum_k = 0
    for i in range(k):
        max_s = 0
        for j in range(s):
            if mat[i][j] > max_s:
                max_s = mat[i][j]
        sum_k += max_s
    return sum_k/total*100

# Get recall
def getRecall(mat, k, s, total, unclassified):
    sum_s = 0
    for i in range(s):
        max_k = 0
        for j in range(k):
            if mat[j][i] > max_k:
                max_k = mat[j][i]
        sum_s += max_k
    return sum_s/(total+unclassified)*100

# Get ARI
def getARI(mat, k, s, N):
    t1 = 0    
    for i in range(k):
        sum_k = 0
        for j in range(s):
            sum_k += mat[i][j]
        t1 += scipy.special.binom(sum_k, 2)
        
    t2 = 0
    for i in range(s):
        sum_s = 0
        for j in range(k):
            sum_s += mat[j][i]
        t2 += scipy.special.binom(sum_s, 2)
        
    t3 = t1*t2/scipy.special.binom(N, 2)
    
    t = 0
    for i in range(k):
        for j in range(s):
            t += scipy.special.binom(mat[i][j], 2)
        
    ari = (t-t3)/((t1+t2)/2-t3)*100
    return ari

# Get F1-score
def getF1(prec, recall):
    return 2*prec*recall/(prec+recall)


# Determine precision, recall, F1-score and ARI for clustering result
#------------------------------------------------------------------

total_clustered = 0

clusters_species = [[0 for x in range(gold_standard_n_clusters)] for y in range(n_clusters)]

for i in range(n_clusters):
    for j in range(gold_standard_n_clusters):
        n = 0
        for k in range(clustered_count):
            if clustered_objects[k] in clusters[i] and clustered_objects[k] in gold_standard_clusters[j]:
                n+=1
                total_clustered += 1
        clusters_species[i][j] = n

results.append("{}{}\n".format("Number of objects available in the clustering result that are present in the gold standard:", total_clustered))
#print("Number of objects available in the clustering result that are present in the gold standard:", total_clustered)

my_precision = getPrecision(clusters_species, n_clusters, gold_standard_n_clusters, total_clustered)
my_recall = getRecall(clusters_species, n_clusters, gold_standard_n_clusters, total_clustered, (gold_standard_count-total_clustered))
my_ari = getARI(clusters_species, n_clusters, gold_standard_n_clusters, total_clustered)
my_f1 = getF1(my_precision, my_recall)

results.append("{}{}\n".format("Precision = ", my_precision))
results.append("{}{}\n".format("Recall = ", my_recall))
results.append("{}{}\n".format("F1-score = ", my_f1))
results.append("{}{}\n".format("ARI = ", my_ari))

#print("\nEvaluation Results:")
#print("Precision =", my_precision)
#print("Recall =", my_recall)
#print("F1-score =", my_f1)
#print("ARI =", my_ari)

#print()
num = clustered_file.split('/')[-1].split('.')[0].split('t')[-1]
file_name = './results_files_batch_folder_2/' + 'result' + num + ".txt"
results_file = open(file_name, "w")
results_file.writelines(results)
results_file.close()