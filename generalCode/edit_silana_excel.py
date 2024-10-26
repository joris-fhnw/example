import numpy as np
from datetime import datetime, date, time


def edit_silana_excel(file):
    file.columns = change_header_name(file.columns)
    file = file.drop(0)  # delete column 0, contains units
    file = file.reset_index(drop=True)  # reset index to 0
    return file

def change_header_name(header):
    """
    header so bearbeiten, dass nur erster Teil verwendet wird
    Bsp.: Temp [°C] wird zu Temp
    so kann direkt mit den headers ohne [] gearbeitet werden

    Parameters
    ----------
    header = array mit allen namen die Verwendet werden sollen
    """
    names = np.empty(len(header), dtype=object)
    z = 0
    for i in header:
        if ' ' in i:
            temp = i.rsplit()
            names[z] = temp[0]
        else:
            names[z] = i  # header

        if "-" in i:
            tmp = i.replace("-", "_")
            names[z] = tmp

        z += 1
    return names

def time_abs_fun(Uhrzeit):
    start_time = datetime(202, 12, 8, Uhrzeit[0].hour, Uhrzeit[0].minute, Uhrzeit[0].second)  # Referenz erstellen, damit mit Time gerechnet werden kann. Mit datetime.time kann nicht gerechnet werden, da kein 0Punkt vorhanden. Dieser muss erstellt werden
    time_abs = np.zeros(len(Uhrzeit))
    for i in range(0, len(time_abs)):
        time_now = datetime(202, 12, 8, Uhrzeit[i].hour, Uhrzeit[i].minute
                            , Uhrzeit[i].second)
        time_abs[i] = np.timedelta64(time_now - start_time, "s").astype(int)
    return time_abs


def delete_time_span(file, header, start, end):
    """
    Werte von start bis ende  und raw-data werden gelöscht


    Parameters
    ----------
    header: Array aller Namen die verwendet werden sollen
    start:  Start der zu löschenden Daten, muss int sein
    end:   Ende der zu löschenden Daten, muss int sein
    """

   # z = 0
    file = file.drop(file.index[start:end])
    file = file.reset_index(drop=True)
    return file

def del_raw(file,col_raw):
    col_raw = file.str.contains("_raw")

    print(col_raw)
    return file

def get_index_of_time(Uhrzeit,year,month,day,hour,minute,second):
    """

    :param Uhrzeit: datetime.datetime, Index wird in diesem Array gesucht
    :param year: 2022
    :param month: 10
    :param day: 1
    :param hour: 17
    :param minute: 10
    :param second: 1
    :return: index
    """
    return np.where(Uhrzeit == np.datetime64(datetime(year, month, day, hour, minute, second)))[0][0]

def get_index_of_time_for_timezone(Uhrzeit,hour,minute,second):
    return Uhrzeit.where((Uhrzeit.dt.hour == hour) & (Uhrzeit.dt.minute == minute) &
                                (Uhrzeit.dt.second == second)).dropna().index[0]