import os
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg") #Qt5Agg


color = {"CO": "green", "CO2":"salmon","NOx":"purple","O2":"royalblue","UHC":"darkorange",
         "P_Kessel_Silana": "darkred"}

def emi_plot(file_name,df,start,end,y1lim,y2lim,save_fig,path):
    """

    :param file_name: titel für den Plot
    :param df: Dataframe aus welchem die Daten geholt werden
    :param start: start Zeitpunkt des pPlots
    :param end: end Zeitpunkt des Plots
    :param save_fig: soll der Plot gespeichert werden?
    :param path: speicherort des Plots
    """
    fig1, ax = plt.subplots()
    plt.title(file_name)
    # Set the format of the x-axis to hh:mm:ss
    xformatter = mdates.DateFormatter('%H:%M:%S')
    #if time_zone != "None":
    #    xformatter.set_tzinfo(timezone(time_zone))
    ax.xaxis.set_major_formatter(xformatter)

    ax.plot(df.Uhrzeit[start:end], df.CO[start:end], color=color["CO"],  label="CO")
    ax.plot(df.Uhrzeit[start:end], df.NO[start:end] , color=color["NOx"], label="NOx")

    ax2 = ax.twinx()
    ax2.plot(df.Uhrzeit[start:end], df.O2[start:end], color=color["O2"], label="O2")


    ax.set_ylabel("[ppm]",)
    ax2.set_ylabel(' [Vol %] und [kW]')

    ax.set_xlabel(' Time [hh,mm,ss]')
    ax.grid(True)
    ax.set_ylim(y1lim)
    ax2.set_ylim(y2lim)

    # fig1.legend(loc='upper right', fontsize=size)

    fig1.tight_layout()

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='best')

    if save_fig:
        fig_formate = "png"
        fig_name = file_name + "_Emissionen" + "." + fig_formate
        fig_path = path + "/" + fig_name
        plt.savefig(fig_path, format=fig_formate)

    plt.show()
