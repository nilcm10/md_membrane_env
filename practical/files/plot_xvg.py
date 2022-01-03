import sys
import pandas as pd
import matplotlib.pyplot as plt

file = str(sys.argv[1])

def read_xvg(file):
    df = pd.DataFrame()
    with open(file) as f:
        for line in f.readlines():
            if '#' not in line and '@' not in line:
                row = [float(x) for x in line.split()]
                df = df.append([row])
    return df

df = read_xvg(file)
df = df.set_index(0, drop=True)

df.plot(subplots=True, grid=True, figsize=(10,5))
fig = plt.gcf()

out = "{}.png".format(file[:-4])
fig.savefig(out, dpi=200)

