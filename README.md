
# Description
## Manhattan Plot by Ankush Yadav
A Manhattan plot is a specific type of scatter plot widely used in genomics to study GWAS results (Genome-Wide Association Study). Each point represents a genetic variant. The X-axis shows its position on a chromosome, the Y-axis tells how much it is associated with a trait.

## Link for downloading the data 
https://www.dropbox.com/s/dxfr1uq20wbdj1d/AUTOMOBILE_SPEEDING_PROPENSITY_GWAS.txt?dl=0

# Import packages
```bash
import pandas as pd
import numpy as np
import seaborn as sns
%matplotlib inline
import adjustText
from adjustText import adjust_text
```
# Upload single nucleotide polymorphism data
```bash
df = pd.read_csv("AUTOMOBILE_SPEEDING_PROPENSITY_GWAS.txt", sep='\t')
df.head()
```
# Convert P-value into -logP
```bash
- np.log(df['Pval'])
```

## Add -logP into the CSV table in order to put on x-axis
```bash
df['-logP'] = - np.log(df['Pval'])
df
```
# Set the position of SNPs on the X-axis
There are certain ways to set the position of SNPs on the x-axis, one method is where just start counting SNPs from position 1 and simply assign a number to each SNP depending on the order that is occured.

Second method used the actual genomic coordinates so have the position of the SNPs but rather than counting from the start of each chromosomes therefore we will group these SNPs according to each chromosome, we group by passing this chromosome column, we get a group by object and then we iterate, we get a tuple .
we use group_df for group dataframe when doing this
```bash
for chrom, group_df in df.groupby('CHR'):
    print(chrom,group_df['POS'].max())
```
Here we ordered the position in chromosomes from biggest to smallest. If we take these maxes then we can add a kind of ruuning total or running position that will start from the zero and we can put some logic here to add on this running position.
```bash
running_pos = 0

for chrom, group_df in df.groupby('CHR'):
    print(chrom, running_pos)
    running_pos += group_df['POS'].max()
```
Here first chromosome start from position 0 relative to itself and then second chromosme start at  249225077 and so on. As we know entire genome is approx 3 billion therefore last chromosome position is at 2828131431. We take these current running position and just add these to the group

```bash
running_pos = 0

for chrom, group_df in df.groupby('CHR'):
    print(chrom, group_df['POS'] + running_pos)
    running_pos += group_df['POS'].max()  
```
These are the adjusted position of all SNPs in respective chromosome. Now we add a list and instead of printing it we will append this running position or series.Then we concatenate this list of all our series one per chromosome and then add all these as new column
```bash
running_pos = 0
cumulative_pos = []

for chrom, group_df in df.groupby('CHR'):
    cumulative_pos.append(group_df['POS'] + running_pos)
    running_pos += group_df['POS'].max()

df['cumulative_pos'] = pd.concat(cumulative_pos)
df
```

## Make an Index as SNP Number and store in CSV table
```bash
df['SNP Number'] = df.index
df
```
# Now visualize the plot using seaborn package
```bash
sns.relplot(data = df, x = 'cumulative_pos', y = '-logP, aspect = 4)
```
Set the color of all the SNPs position according to their position on the chromosome.
```bash
sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1')
```
If we want, we can also use SNP Number on the x-axis
```bash
sns.relplot(data = df.sample(100000) , x = 'SNP Number', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1')
```
We will change the linewidth, size , and turnoff the legend and store the result as 'g'
```bash
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1', linewidth = 0, s = 6, legend = None)
g
```
We will label x-axis as chromsome
```bash
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1', linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
```
Now we want each chromosome to be labeled on the x-axis. we can easily identify the position of each chromosome on the x-axis by taking data frame and group it by chromosome and take cumulative position and then take median
```bash
df.groupby('CHR')['cumulative_pos'].median()
```
We will add this position on the x-axis by using tick.label code.
```bash
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1', linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
```
Another method of  tick labeling
```bash
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = 'Set1', linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
```
if you want other colors such as green and yellow. we have 22 chromosome, so we want 11 green and 11 yellow at alternate location.
```bash
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = ['green', 'yellow']*11, linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
```
Alternative method to set two color in at alternate position. if you want to give odd chromosome- a value A and even chromosomes  value say B. 
```bash
df['color group'] = df['CHR'].apply(lambda x : 'A' if x % 2 == 0 else 'B')
df

g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow'], linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
```
If we want to highlight only SNPs on chromosome 4 with specific color and denoted it as C.
```bash

df['color group'] = df['CHR'].apply(lambda x : 'A' if x % 2 == 0 else 'B')
df.loc[df['CHR'] == 4, 'color group'] = 'C'

g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow', 'red'], linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
```
if we want to highlight only specific position.
```bash
df['color group'] = df['CHR'].apply(lambda x : 'A' if x % 2 == 0 else 'B')
df.loc[(df['CHR'] == 3) & (df['POS'] < 50000000, 'color group'] = 'C'

g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow', 'red'], linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].medium())
g.ax.set_xticklables(df['CHR'].unique())
```
If we want to set cutoff for significance .
```bash

df['color group'] = df['CHR'].apply(lambda x : 'A' if x % 2 == 0 else 'B')
df.loc[(df['CHR'] == 3) & (df['POS'] < 50000000, 'color group'] = 'C'

g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow', 'red'], linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].medium())
g.ax.set_xticklables(df['CHR'].unique())
g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow', 'red'], linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
g.ax.axhline(12, linestyle='--', linewidth = 1)
```
 We can also plot manhattan plot for each bases.
 ```bash

 df['color group'] = df['CHR'].apply(lambda x : 'A' if x % 2 == 0 else 'B')
df.loc[(df['CHR'] == 3) & (df['POS'] < 50000000, 'color group'] = 'C'

g = sns.relplot(data = df.sample(100000) , x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'color group', palette = ['green', 'yellow', 'red'], linewidth = 0, s = 6, legend = None, row = 'A1')
for ax in g.axes.flat:
    ax.set_label('Chromosome')
    ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
    ax.set_xticklables(df['CHR'].unique())
    ax.axhline(12, linestyle='--', linewidth = 1)

```

If we want to annotate particular set of markers or position. set some threshold say -logP=20.
```bash
my_data[my_data['-logP'] > 20]
```
Here we use apply method  where we use lambda expression and then annotate method and p is used an argument
```bash
my_data = df.sample(100000)
g = sns.relplot(data = my_data, x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = ['green', 'yellow']*11, linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
my_data[my_data['-logP'] > 20].apply(lambda p : g.ax.annotate(p['MarkerName'], (p['cumulative_pos'], p['-logP'])), axis=1)

```
To adjust the annotation.

```bash
import adjustText
from adjustText import adjust_text

my_data = df.sample(100000)
g = sns.relplot(data = my_data, x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = ['green', 'yellow']*11, linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
annotations = my_data[my_data['-logP'] > 20].apply(lambda p : g.ax.annotate(p['MarkerName'], (p['cumulative_pos'], p['-logP'])), axis=1).to_list()
adjust_text(annotations)
```
if we want to clarify more so that we can easily analyze and see annotations. 
```bash
import adjustText
from adjustText import adjust_text

my_data = df.sample(100000)
g = sns.relplot(data = my_data, x = 'cumulative_pos', y = '-logP', aspect = 4, hue = 'CHR', palette = ['green', 'yellow']*11, linewidth = 0, s = 6, legend = None)
g.ax.set_label('Chromosome')
g.ax.set_xticks(df.groupby('CHR')['cumulative_pos'].median())
g.ax.set_xticklables(df['CHR'].unique())
annotations = my_data[my_data['-logP'] > 20].apply(lambda p : g.ax.annotate(p['MarkerName'], (p['cumulative_pos'], p['-logP'])), axis=1).to_list()
adjust_text(annotations, arrowprops = {'arrowstyle' : '>', 'color' : 'purple'})
```











