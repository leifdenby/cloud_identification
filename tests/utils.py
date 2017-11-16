import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plot

def plot_labels(x, y, labels, filename):
    plot.figure()

    # plot.pcolormesh(x, y, labels)
    plot.imshow(labels, interpolation='nearest', extent=[x.min(), x.max(), y.min(), y.max()])
    plot.gca().set_aspect(1)
    plot.colorbar()

    plot.savefig(filename)
