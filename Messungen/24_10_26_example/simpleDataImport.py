import pandas as pd
import os
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg") #Qt5Agg


"Daten Import"
df = pd.read_csv("exampleData.txt",delimiter="\t",encoding = 'unicode_escape',
                 header=1,low_memory=False)
df = df.drop(0) # löscht erste Zeile in txt
df = df.reset_index(drop=True)
df.Uhrzeit = pd.to_datetime(df.Uhrzeit,dayfirst=True) #Uhrzteit to datetime
df = df.astype({col: float for col in df.columns[1:]}) #Datentyp in float ändern (alle ausser Uhrzeit)

"Plot"
fig1, ax = plt.subplots()
# Set the format of the x-axis to hh:mm:ss
xformatter = mdates.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(xformatter)

ax.plot(df.Uhrzeit, df.CO, color="red",  label="CO")
ax2 = ax.twinx()
ax2.plot(df.Uhrzeit, df.O2, color="blue", label="O2")
ax.set_ylabel("[ppm]",)
ax2.set_ylabel(' [Vol %] und [kW]')
ax.set_xlabel(' Time [hh,mm,ss]')
ax.grid(True)
ax.set_ylim([0,100])
plt.show()
