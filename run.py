from multiprocessing import Process
import subprocess as sp
import numpy as np

### !! NOT UP TO DATE !! ###


# clean all preCICE logs
sp.call(['sh', './clean.sh'])


# def all functions, one per type of participant 
def cell(i):
    path = "cell_" + str(i) + "/forced_convection_" + str(i) + ".py"
    exec(open(path).read())


# main 
if __name__ == "__main__":

    # start all participants

    #cells 
    num_cells = 2
    cells = None 
    i = 1
    while i <= num_cells:
        cells = np.append(cells, Process(target=cell, args=(i,)))
        cells[i].start()
        i += 1



    # join all participants 

    #cells 
    i = 1
    while i <= num_cells:
        cells[i].join()
        i += 1

