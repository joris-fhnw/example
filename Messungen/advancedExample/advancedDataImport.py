import pandas as pd
import os
from generalCode import plots
from generalCode import edit_silana_excel as ed
from generalCode import combsutionFunctions as cf
import matplotlib
matplotlib.use("Qt5Agg") #Qt5Agg // TkAgg



show_plot = True
save_fig=True


" Daten Import"
df = pd.read_csv("exampleData.txt",delimiter="\t",encoding = 'unicode_escape',
                 header=1,low_memory=False)
df = df.drop(0)
df = df.reset_index(drop=True)
df.Uhrzeit = pd.to_datetime(df.Uhrzeit,dayfirst=True)
df = df.astype({col: float for col in df.columns[1:]})

"Mittelwertbildung"
#Datum des Messtages (Datum im txt)
year = 2023
month = 7
day = 31

start = ed.get_index_of_time(df.Uhrzeit, year, month, day, 15, 50, 0)
stopp = ed.get_index_of_time(df.Uhrzeit, year, month, day, 17, 35, 0)

mean_values = df[start:stopp].mean()
mean_values = mean_values.round()

"Berechnungen"
# Wood properties
w = 0.08  # [kg wasser/ kg holz nass]
yc, yh, yo, yn, ys = 0.51, 0.06, 0.43, 0, 0
gamma = ["C", yc, "H", yh, "O", yo, "N", yn, "S", ys]
pellet = cf.Wood("Pellet", gamma, w)  # in pellet are all the needed combustion calculation stored
lmin = round(pellet.lmin(),2)
lam = 2  #lambda der Verbrennung
o2_calc = round(pellet.O2(lam)*100,1)

print(f"Mindestlufmtenge für ein Pellet: {lmin} kg Luft/ kg Pellet\n"
      f"bei einem Lamba von {lam} beträgt der O2-Gehalt im Abgas: {o2_calc} VOl%")

"Plot"
if show_plot:
    # Neuen Pfad erstellen --> Speicherort des Plots
    path = os.getcwd()
    path = path + "/" +'All'+'_figures'
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

    file_name = os.path.basename(path)

    start=0
    end=df.O2.size
    y1lim = [0,300]
    y2lim=[0,21]
    #Aufrufen der Plot Funktion
    plots.emi_plot(file_name,df,start,end,y1lim,y2lim,save_fig,path)
