import sys
import pandas as pd
import matplotlib.pyplot as plt

def from_letter_to_value(letter):
    return letter2value[letter]

file = str(sys.argv[1])

rows = []
with open(file, "r") as f:
    for i, line in enumerate(f.readlines()):
        rows.append([*line[:-1]])
df = pd.DataFrame(rows)

letters = list(set(df.values.flatten()))
values = [float(x) for x in range(len(letters))]
letter2value = dict(zip(letters, values))
for i in range(len(df)):
    df.loc[i] = df.loc[i].apply(from_letter_to_value)

fig, ax = plt.subplots(figsize=(20,15))
ax.imshow(df.values.astype(float))
out = "{}.png".format(file[:-4])
fig.savefig(out, dpi=400)
