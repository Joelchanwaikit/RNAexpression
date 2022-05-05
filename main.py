import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

def RNAexpression_analyser():
    # 1 Read the file. User to define path
    df_unclean = pd.read_csv(r"C:\Users\joelc\OneDrive\Desktop\Input\RT-qPCR_input.csv")

    # 2 Load in Data
    df = df_unclean[df_unclean["Content"] == "Unkn"]
    col_sampletype = list(df.Target + ' ' + df.Sample)
    for z in range(len(col_sampletype)):
        col_sampletype[z] = col_sampletype[z].lower()
    col_cqvalues = list(df.Cq)

    # 3 Sort by sample type and the cq values
    Clean_data = {}
    for i in range(len(col_cqvalues)):
        if col_sampletype[i] in Clean_data.keys():
            Clean_data[col_sampletype[i]] = Clean_data[col_sampletype[i]] + [col_cqvalues[i]]
        else:
            Clean_data[col_sampletype[i]] = [col_cqvalues[i]]

    # 4 Have the sample type with the Cq mean and Std
    Analysed_data = {}
    for key in Clean_data:
        Values = np.array(Clean_data[key])
        Cq_mean = np.mean(Values)
        Cq_std = np.std(Values)
        Analysed_data[key] = [Cq_mean, Cq_std]
    # Export Mean Cq Data
    sample = []
    Cq = []
    Std = []
    for keys in Analysed_data:
        sample.append(keys)
        Cq.append(Analysed_data[keys][0])
        Std.append(Analysed_data[keys][1])
    data = {"Samples": sample, "CQ": Cq, "Std": Std}

    # 5 Ask what is the housekeeping gene. sort in to housekeeping and target
    Housekeeping = input("Housekeeping Gene? Maximum 1").strip().lower()
    Housekeeping_no = len(Housekeeping)

    Housekeeping_details = {}
    Target_details = {}
    for key in Analysed_data:

        if Housekeeping in key:
            Housekeeping_details[key.replace(Housekeeping, '')] = Analysed_data[key]
        else:
            Target_details[key] = Analysed_data[key]

    # 6 Ask for target gene names, seperated by spaces
    TargetsID = list(input("Name of target genes, seperated by comma").lower().split(','))

    # 7 Ask for sample  names, seperated by spaces
    SamplesID = list(input("Name of Samples , seperated by comma").lower().split(','))

    # 8 Delta Cq and Std deviation
    deltacq_table = {}
    for keys in Target_details:
        Target_cq_mean = Target_details[keys][0]
        for i in range(len(SamplesID)):
            if SamplesID[i] in keys:
                Sample = " " + SamplesID[i]
        Housekeeping_mean = Housekeeping_details[Sample][0]
        deltacq = Target_cq_mean - Housekeeping_mean
        std = Target_details[keys][1] + Housekeeping_details[Sample][1]
        deltacq_table[keys] = [deltacq, std]

    # 9 DeltaDelta Comparison
    sample = []
    deltaCq = []
    Std = []
    for keys in deltacq_table:
        sample.append(keys)
        deltaCq.append(deltacq_table[keys][0])
        Std.append(deltacq_table[keys][1])

    arbitary_value = max(deltaCq)
    deltadeltacq = []
    relative_expression = []
    err_1 = []
    err_2 = []
    for i in range(len(sample)):
        deltadeltacq.append(deltaCq[i] - arbitary_value)
        relative_expression.append(math.pow(2, arbitary_value - deltaCq[i]))
        err_1.append(math.pow(2, arbitary_value - deltaCq[i] - Std[i]))
        err_2.append(math.pow(2, arbitary_value - deltaCq[i] + Std[i]))

    data = {"Samples": sample, "DeltaCQ": deltaCq, "Std": Std, "DDCQ": deltadeltacq,
            "Relative_Expression": relative_expression, "err_1": err_1, 'err_2': err_2}
    dfexport = pd.DataFrame(data)
    dfexport.to_csv('validation.csv')


def RNAplotter():
    # 1 Read files from Analysis Module
    df = pd.read_csv('validation.csv')

    # 2 Load Data
    sampleID = list(df.Samples)
    sampleID = [i.split(' ') for i in sampleID]
    DDCQ = list(df.Relative_Expression)
    lower_err = list(df.err_1)
    upper_err = list(df.err_2)

    # 3 Get info on what to plot
    samples = list(input("Name of Samples , seperated by comma").lower().split(','))
    targets = list(input("Name of target genes, seperated by comma").lower().split(','))

    # 4 Get error bars and plot each bar
    upper1 = []
    lower1 = []
    w = 0.2
    for i in range(len(DDCQ)):
        lower1.append(upper_err[i] - DDCQ[i])
        upper1.append(DDCQ[i] - lower_err[i])

    for i in range(len(samples)):
        to_plot = []
        for j in range(len(sampleID)):
            if samples[i] == sampleID[j][1]:
                label = samples[i]
                to_plot.append(int(j))
        ddcq_toplot = []
        targets_toplot = []
        upper_err_toplot = []
        lower_err_toplot = []
        for k in range(len(to_plot)):
            ddcq_toplot.append(DDCQ[to_plot[k]])
            targets_toplot.append(sampleID[to_plot[k]][0])
            upper_err_toplot.append(upper1[to_plot[k]])
            lower_err_toplot.append(lower1[to_plot[k]])
        asymmetric_error = np.array(list(zip(upper_err_toplot, lower_err_toplot))).T
        x = list(np.arange(len(targets_toplot)))
        x = [x + w * i for x in x]

        plt.bar(x, ddcq_toplot, width=w, label=F"{label}")
        plt.errorbar(x=x, y=ddcq_toplot, yerr=asymmetric_error, c='k', fmt=' ', elinewidth=1, capsize=3)

    # 5 Define the graph items

    plt.xticks(np.array(x) - w * (1.5), targets_toplot)
    plt.legend()
    plt.title(F"Relative Gene Expression")
    plt.xlabel("Genes")
    plt.ylabel("Relative Expression")

    plt.show()

    repeat = input("Another Graph? (y/n)")
    if repeat == 'y':
        RNAplotter()


if __name__ == '__main__':
    print("a: Analyser , b: Plotter")
    program = input("Which module would you like to use? (a/b)")
    if program == 'a':
        RNAexpression_analyser()
        RNAplotter()
    elif program == 'b':
        RNAplotter()
    else:
        print("Invalid input")

