import numpy as np
import matplotlib.pyplot as plt

q = 2
r = 5
v = 8
o = 11
l = 14
time = 17

colors = dict({1:"#770000", 2: "#992222", 3:"#aa4444", 4:"#cc6666", 5:"#ff8888", 256: "#aa22aa", 16: "#22aaaa"})

matrix = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

file = open("timings_extension.dat", "r")
for line in file:
    words = line.split()
    matrix = np.vstack([matrix, [int(words[q]), int(words[r]), int(words[v]), int(words[o]), int(words[l]), np.log(0.00001 + float(words[time]))/np.log(2.0)]])

# get all q's
rr = list(set([matrix[i,1] for i in range(0, np.size(matrix,0))]))

for r in rr:
    if r == 0:
        continue

    indices = [i for i in range(0, np.size(matrix,0)) if matrix[i,1] == r]
    mat = matrix[indices,:]

    # get all o's for this r
    oo = list(set(mat[i,3] for i in range(0, np.size(mat,0))))
    oo.sort(key=int)
    for o in oo:
        # do we want to plot this line? maybe avoid clutter
        if r == 2:
            continue

        indices = [i for i in range(0, np.size(mat,0)) if mat[i,3] == o]
        M = mat[indices,:]

        # get all l's for this pair (o, q)
        ll = [M[i,4] for i in range(0, np.size(M, 0))]
        ll = list(set(ll))
        ll.sort(key=int)
        if len(ll) == 0:
            continue

        # for all l, compute the average timing and add it to the list
        avgs = [0.0] * len(ll)
        for j in range(0, len(ll)):
            ar = [M[i,5] for i in range(0, np.size(M, 0)) if M[i,4] == ll[j]]
            avgs[j] = np.mean(ar)

        # add line to plot
        plt.plot(ll, avgs, 'o-', color=colors[int(r)], linewidth=2, markersize=2)
        location = (ll[0]-0.2, avgs[0])
        if r == 1 and o in range(8,16):
            location = (ll[0]-0.2, avgs[0]-0.5)
        plt.annotate('$r = %i, m = %i$' % (int(r), int(o)), xy=location, horizontalalignment='right', verticalalignment='center', color=colors[int(r)], fontsize=13, fontweight='bold')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel('degree of circulancy $\ell$')
plt.ylabel('logarithm base 2 of solving time (seconds)');
plt.xlim(-2.5, 13)
plt.ylim(-10, 15)
plt.show()

#from matplotlib2tikz import save as tikz_save
#tikz_save('plot.tex') # don't forget to add '/pgf/number format/1001 sep={}' in the [] prefix of the axis environment

