import numpy as np
import math
import random
import copy
import statistics as stat
import matplotlib.pyplot as plt
from matplotlib import cm


### Question 1 - write function for K means algorithm ###
def Kmeans(X, start, nIter):

    # Input for the function:
    # X - a matrix of the points. Each row is a signle point. Each column is a dimension.
    # start - a matrix of the centers. Each row is a signle center. Each column is a dimension. K = num of rows.
    # nIter - number of iterations for the algorithm.

    nPoints = len(X[0,:])
    K = len(start[0,:])
    centers = copy.copy(start)
    cluster = np.zeros((nPoints,), dtype = int)

    # loop to create iterations
    for iteration in range(0, nIter):

        # create new figure
        plt.figure()
        
        # Run on all points
        for pointNum in range(0, nPoints):

            # Reset distance vector for point
            disVec = np.zeros((K,), dtype = int)

            # Run on all K centers
            for i in range(0, K):
                
                # Calculate distance for current center
                disVec[i] = math.sqrt((X[0, pointNum] - centers[0, i]) ** 2 + (X[1, pointNum] - centers[1, i]) ** 2)

            # Reassign point to center based on minimal distance
            cluster[pointNum] = np.argmin(disVec)

        # Compute new centers and plot scatters for new clusters and centers
        colors = iter(cm.rainbow(np.linspace(0, 1, K)))
        
        for i in range(0,K):

            x = X[0, cluster == i]
            y = X[1, cluster == i]

            # calculate new center
            centers[0, i] = stat.mean(x)
            centers[1, i] = stat.mean(y)

            # plot scatter for new center
            plt.scatter(x, y, color=next(colors), s = 40)
            plt.scatter(centers[0, i], centers[1, i], color = 'black', s = 15)
            title = "Iteration No. %i - kMeans = %d" % (iteration+1, K)
            plt.title(title)
            plt.xlabel("x axis")
            plt.ylabel("y axis")

    # On last figure plot initial centers
    plt.scatter(start[0, :], start[1, :], color = 'black', s = 20, marker = '*')

    return[centers, cluster]



### Question 2 ###

X = np.array([[2, 2, 8, 5, 7, 6, 1, 4], [10, 5, 4, 8, 5, 4, 2, 9]])
start = np.array([[2, 5, 1],[10, 8, 2]])
nIter = 3

[centers, cluster] = Kmeans(X,start,nIter)

print(centers)
print(cluster)




### Question 3 - run for 3 normal distributions as examples ###

Kvec = [2,3,4]
nIter = 10

# Initialize values for normal distribution
mu = [1, 11, 22]
sigma = 3
normDisX = []
normDisY = []

# Create three normal distributions
for i in range(0, 3):

    disX = np.random.normal(mu[i], sigma, 100)
    disY = np.random.normal(mu[i], sigma, 100)

    normDisX = np.append(normDisX, disX)
    normDisY = np.append(normDisY, disY)

# Add x axis and y axis to matrix of normal distribution
normDis = np.vstack((normDisX, normDisY))

# Run algorithm for each K value
for K in Kvec:

    xCenter = []
    yCenter = []

    # Choose K random centers
    for k in range(0, K):

        randX = random.randint(0,20)
        randY = random.randint(0,20)
        xCenter = np.append(xCenter, randX)
        yCenter = np.append(yCenter, randY)

    initCenter = np.vstack((xCenter,yCenter))

    # Run k means algorithm
    [centers, cluster] = Kmeans(normDis, initCenter, nIter)
