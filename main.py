import os

path = './ground_thruth_files_batch_folder_2'
files = os.listdir(path)
folder = []
for f in files:
    folder.append(f)

python_cmd = "python evaluate.py "
clustered = "--clustered ./clustering_results_batch_folder_2/clustering_result"
goldstandard = " --goldstandard ./ground_thruth_files_batch_folder_2/ground_thruth"


for i in range(1,len(folder)):
    command = "{}{}{}{}{}{}{}".format(python_cmd, clustered, i, ".csv", goldstandard, i, ".csv" )
    os.system(command)

