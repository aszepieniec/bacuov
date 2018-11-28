import numpy as np
import matplotlib.pyplot as plt

q = 2
v = 5
o = 8
l = 11
time = 14

colors = dict({2:"#aa2222", 7:"#22aa22", 31:"#22aaaa", 251:"#aaaa22", 256: "#aa22aa", 16: "#22aaaa"})

matrix = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])

file = open("timings.dat", "r")
for line in file:
    words = line.split()
    matrix = np.vstack([matrix, [int(words[q]), int(words[v]), int(words[o]), int(words[l]), np.log(0.00001 + float(words[time]))/np.log(2.0)]])

# get all q's
qq = list(set([matrix[i,0] for i in range(0, np.size(matrix,0))]))

for q in qq:
    if q == 0:
        continue

    # which plot to make? set == or !=
    if int(q) % 2 == 0:
        continue

    indices = [i for i in range(0, np.size(matrix,0)) if matrix[i,0] == q]
    mat = matrix[indices,:]

    # get all o's for this q
    oo = list(set(mat[i,2] for i in range(0, np.size(mat,0))))
    oo.sort(key=int)
    for o in oo:
        # do we want to plot this line? maybe avoid clutter
        if q == 16 and o == 20:
            continue
        elif q == 16 and o == 12:
            continue
        elif q == 16 and o == 6:
            continue
        elif q == 16 and o == 9:
            continue
        elif q == 2 and o == 21:
            continue
        elif q == 7 and o == 8:
            continue
        elif q == 256 and o == 20:
            continue
        elif q == 2 and o == 20:
            continue
        elif q == 256 and o == 15:
            continue
        elif q == 256 and o == 8:
            continue
        elif q == 256 and o == 4:
            continue
        elif q == 256 and o == 25:
            continue

        indices = [i for i in range(0, np.size(mat,0)) if mat[i,2] == o]
        M = mat[indices,:]

        # get all l's for this pair (o, q)
        ll = [M[i,3] for i in range(0, np.size(M, 0))]
        ll = list(set(ll))
        ll.sort(key=int)
        if len(ll) == 0:
            continue

        # for all l, compute the average timing and add it to the list
        avgs = [0.0] * len(ll)
        for j in range(0, len(ll)):
            ar = [M[i,4] for i in range(0, np.size(M, 0)) if M[i,3] == ll[j]]
            avgs[j] = np.mean(ar)

        # add line to plot
        plt.plot(ll, avgs, 'o-', color=colors[int(q)], linewidth=2)
        location = (ll[0]-0.2, avgs[0])
        if q == 7 and o == 8:
            location = (ll[0]-0.2, avgs[0]-0.5)
        elif q == 251 and o == 3:
            location = (ll[0]-0.2, avgs[0]-0.4)
        elif q == 16 and o == 16:
            location = (ll[0]-0.2, avgs[0]-0.4)
        elif q == 16 and o == 15:
            location = (ll[0]-0.2, avgs[0]-0.4)
        elif q == 256 and o == 12:
            location = (ll[0]-0.2, avgs[0]-0.6)
        elif q == 256 and o == 27:
            location = (ll[0]-0.2, avgs[0]+0.4)
        plt.annotate('$q = %i, m = %i$' % (int(q), int(o)), xy=location, horizontalalignment='right', verticalalignment='center', color=colors[int(q)], fontsize=13, fontweight='bold')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel('degree of circulancy $\ell$')
plt.ylabel('logarithm base 2 of solving time (seconds)');
plt.xlim(-2.5, 13)
plt.ylim(-10, 15)
plt.show()

#from matplotlib2tikz import save as tikz_save
#tikz_save('plot.tex') # don't forget to add '/pgf/number format/1001 sep={}' in the [] prefix of the axis environment

