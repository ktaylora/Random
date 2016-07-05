import pandas as pd
import matplotlib.pyplot as plt

t = pd.read_csv("/Users/ktaylora/Downloads/2013/checklists_50_record_sample.csv")
# t.iloc[:, 2].plot(kind='hist')  # LATITUDE
plt.figure()
pd.DataFrame.hist(t, column="LONGITUDE")
