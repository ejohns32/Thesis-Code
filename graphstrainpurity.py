#! /usr/bin/env python3

import os
import sys
import matplotlib
from matplotlib import pyplot
import numpy

import pandas
from collections import Counter

import locale
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


gold='#FADA5E'
green='#0A7951'
green_gold ='#82AA58'
color_vals = ['purple', 'orange', 'red', green, 'blue', 'black', 'magenta', 'red']
color_vals = ['purple', 'red', 'green', 'yellow',]
marker_vals = ['o','v','p','H','D','s','8','<','>']
marker_vals = ['o','v','s','D']
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'legend.loc' : 'best',
          'lines.color' : green
          }
matplotlib.rcParams.update(params);
default_color = green
def print_counts_table(cluster_counts, caption="", reference="", file=sys.stdout):
    print("\\begin{table}[ht!]", file=file)
    print("\\centering", file=file)
    print("\\begin{tabular}{|c|c|c|}", file=file)
    print("\\hline", file=file)
    print("\\bf %s &\\bf %s &\\bf %s\\\\ \\hline \\hline" % ("\\Spec{}","Count", "Proportion"), file=file)
    sorted_counts = [(k,v) for k,v in cluster_counts.items()]
    sorted_counts.sort(reverse=True, key=lambda x: x[1])
    total = sum([v for k,v in cluster_counts.items()])
    for species, count in sorted_counts:
        print("%s & %s & %.3f\\\\ \\hline" % (species, count, count / total), file=file)
    print("\\hline", file=file)
    print("\\bf %s & \\bf %s & \\bf %.3f\\\\ \\hline" % ("Total", total, 1.000), file=file)
    print("\\end{tabular}", file=file)
    print("\\label{%s}" % reference, file=file)
    print("\\caption{%s}" % caption, file=file)
    print("\\end{table}", file=file)
    print(file=file)
    return
def table_purities(min_pts, cluster_counts):
    sorted_clusters = sorted(cluster_counts, key = lambda x: sum([v for k,v in x.items()]))
    total_clusters = len(cluster_counts)
    i = 0
    with open("cluster-counts-%d.tex" % min_pts, "w") as file:
        for cluster in cluster_counts:
            print("building table %d of %d" % (i, total_clusters))
            i += 1
            if i % 10 == 0:
                print("\\clearpage", file=file)
            print("%%%%% CLUSTER {} of {} NEIGHBORS = {} %%%%%".format(i,total_clusters,min_pts), file=file)
            caption = "Cluster %d of %d for \\minneigh{}=%d" % (i, total_clusters, min_pts)
            reference = "tab:cluster:%d:%d" % (i, min_pts)
            print_counts_table(cluster, reference=reference, caption=caption, file=file)
    return
def table_species(species_counts, filename = "species_table.tex"):
    begin_table = "\\begin{tabular}{c|c|c}\n"
    header = "\\bf %s &\\bf %s \\\\ \\hline \\hline" % ("\\Spec{}", "Number of \\Isols{}")
    end_table = "\\end{tabular}"
    entries = []
    for k,v in species_counts.items():
        entries.append((k,int(v)))
    string = ""
    string += begin_table
    string += header
    entries.sort(key=lambda x : x[1], reverse=True)
    for k,v in entries:
        string += "%s & %s \\\\ \\hline\n" % (k,v)
    string += "\n"
    string += end_table
    print(string)
    with open(filename, "w") as file:
        print(string, file=file)
    return
def plot_hist_species(species_counts, title='Number of Isolates for Each Species ({} Total Isolates)', filename='species_hist.pdf'):
    species_name = []
    species_num = []
    species_list = []
    for k,v in species_counts.items():
        species_name.append(k)
        species_num.append(v)
        species_list.extend([k * v])
    sorted_species_name = []
    sorted_species_num = []
    sorted_species_list = []
    for s,c in sorted(zip(species_name, species_num), key = lambda x: x[1], reverse=True):
        sorted_species_name.append(s)
        sorted_species_num.append(c)
        sorted_species_list.extend([s * c])
    df = pandas.DataFrame({'Species': pandas.Series(sorted_species_name), 'Count':pandas.Series(sorted_species_num) })
    #df = pandas.DataFrame.from_dict(species_counts, orient='index')
    #print('DataFrame               :{}'.format(df))
    #print('DataFrame.columns       :{}'.format(df.columns))
    #print('DataFrame.values        :{}'.format(df.values))
    #print('DataFrame.columns.values:{}'.format(df.columns.values))
    #print('DataFrame.Axes          :{}'.format(df.axes))

    ax = df.plot(kind='bar', title=title.format(sum(species_num)), width=.8, linewidth=.5, color=green, legend=False)
    figure = ax.get_figure()
    #ax = df.apply(pandas.value_counts).plot(kind='bar', subplots=True, title=title, width=.8, linewidth=.5, color=green, legend=False)
    pyplot.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off')         # ticks along the top edge are off
    for p in ax.patches:
        #ax.annotate(locale.format("%d",p.get_height(), grouping=True), (p.get_x(), p.get_height() + 10), rotation=90, verticalalignment='bottom', horizontalalignment='left', fontsize=8)
        ax.annotate('%d'%p.get_height(), (p.get_x(), p.get_height() + 10), rotation=90, verticalalignment='bottom', horizontalalignment='left', fontsize=8)
    pyplot.xticks(range(len(sorted_species_name)), sorted_species_name, fontsize=8)
    pyplot.yticks(numpy.arange(0, 1901, 100))
    ax.get_yaxis().set_tick_params(direction='out')
    ax.set_ylim(0, 1901)
    ax.set_xlabel("Species")
    ax.set_ylabel("Number of Isolates")
    figure.savefig(filename, bbox_inches='tight')
    pyplot.close()
    return None

def plot_two_scatter(x, y1, y2, title, x_label, y_label, filename='scatter.png',
        x_low=None, x_high=None, x_ticks=None, y_low=None, y_high=None, y_ticks=None, marker_size=20, alpha=1, y1_label='y1_label', y2_label='y2_label'):
    print('## SCATTER: {}'.format(filename))
    print('x : {}'.format(x))
    print('y1: {}'.format(y1))
    print('y2: {}'.format(y2))
    axes = pyplot.gca()

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    if x_ticks:
        pyplot.xticks(numpy.arange(x_low, x_high, x_ticks))

    if y_low < y_high:
        axes.set_ylim(y_low, y_high)
    else:
        axes.set_ylim(y_high, y_low)
    if y_ticks:
        pyplot.yticks(numpy.arange(y_low, y_high, y_ticks))
    axes.set_ylim(y_low, y_high)
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.scatter(x, y1, c=green, linewidths=.5, s=marker_size, alpha=alpha, label=y1_label)
    pyplot.scatter(x, y2, c=gold,  linewidths=.5, s=marker_size, alpha=alpha, label=y2_label)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Legend', ncol=2)
    pyplot.savefig(filename, bbox_inches='tight');
    pyplot.close()

def plot_scatter(x, y, title, x_label, y_label, filename='scatter.png',
        x_low=None, x_high=None, x_ticks=None, y_low=None, y_high=None, y_ticks=None, marker_size=20, alpha=1):
    print('## SCATTER: {}'.format(filename))
    axes = pyplot.gca()

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    if x_ticks:
        pyplot.xticks(numpy.arange(x_low, x_high, x_ticks))

    if y_low is None:
        y_low = min(y[0],y[-1]) - 1
    if y_high is None:
        y_high = max(y[0],y[-1]) + 1
    if y_low < y_high:
        axes.set_ylim(y_low, y_high)
    else:
        axes.set_ylim(y_high, y_low)
    if y_ticks:
        pyplot.yticks(numpy.arange(y_low, y_high, y_ticks))
    axes.set_ylim(y_low, y_high)
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.scatter(x, y, c=default_color, linewidths=.5, s=marker_size, alpha=alpha)
    pyplot.title(title)
    pyplot.savefig(filename, bbox_inches='tight');
    pyplot.close()

def plot_lower_upper_unpure_scatter_filled_stack(x, lower, mid_low, mid_high, total, title,
        filename='cluster_stacked_unpure.pdf', x_low=None, x_high=None, y_low=None,
        x_label='Minimum Neighbors', y_label='Number of Isolates',
        lower_label='Misses', mid_low_label='Hits', mid_high_label='Pure Points',
        total_label='Noise',
        legend_title='Isolate Clustering', flip_colors=True,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    print('lower   :{}'.format(lower))
    print('mid-low :{}'.format(mid_low))
    print('mid-high:{}'.format(mid_high))
    print('total   :{}'.format(total))

    assert(all(l + ml + mh <= t for l,ml,mh,t in zip(lower, mid_low, mid_high, total)))
    # COVERAGE VALUES
    with open("coverage-values.tex", "w") as file:
        print("\\begin{table}[ht!]", file=file)
        print("\\setlength{\\arrayrulewidth}{1pt}", file=file)
        print("\\newcolumntype{n}{>{\columncolor{cyan}}r}", file=file)
        print("\\newcolumntype{p}{>{\columncolor{cpgreen}}>{\color{white}}r}", file=file)
        print("\\newcolumntype{h}{>{\columncolor{cpgreengold}}r}", file=file)
        print("\\newcolumntype{m}{>{\columncolor{cpgold}}r}", file=file)
        print("\\centering", file=file)
        print("\\begin{tabular}{|r|m|h|p|n|r|}", file=file)
        print("\\hline", file=file)
        print("\\bf %s & \\bf %s & \\bf %s & \\bf %s & \\bf %s & \\bf %s" %
                ("\\minneigh{}",
                    lower_label,
                    mid_low_label,
                    mid_high_label,
                    total_label,
                    "Total"),
                    file=file)
        print("\\\\", file=file)
        print("\\hline", file=file)
        print("\\hline", file=file)
        for (i,l,ml,mh,h) in zip(x,lower, mid_low, mid_high, total):
            clustered = sum((l,ml,mh))
            print("\\bf %d & %d & %d & %d & %d & %d" % (
                i,
                l,
                ml,
                mh,
                h - clustered,
                h),
            file=file)
            print("\\\\", file=file)
            print("\\hline", file=file)
        print("\\end{tabular}", file=file)
        print("\\caption{The number of \\isols{} that ended up in clusters where their \\spec{} is the minority (%s), the dominant (%s), and make up the entirety of the cluster (%s) or were categorized as %s.}" % 
                (lower_label,mid_low_label,mid_high_label, total_label), file=file)
        print("\\label{tab:coverage:values}", file=file)
        print("\\end{table}", file=file)
    # COVERAGE PERCENT OF TOTAL
    with open("coverage-percent-total.tex", "w") as file:
        print("\\begin{table}[ht!]", file=file)
        print("\\setlength{\\arrayrulewidth}{1pt}", file=file)
        print("\\newcolumntype{n}{>{\columncolor{cyan}}r}", file=file)
        print("\\newcolumntype{p}{>{\columncolor{cpgreen}}>{\color{white}}r}", file=file)
        print("\\newcolumntype{h}{>{\columncolor{cpgreengold}}r}", file=file)
        print("\\newcolumntype{m}{>{\columncolor{cpgold}}r}", file=file)
        print("\\centering", file=file)
        print("\\begin{tabular}{|r|m|h|p|n|r|}", file=file)
        print("\\hline", file=file)
        print("\\bf %s & \\bf %s & \\bf %s & \\bf %s & \\bf %s & \\bf %s" %
                ("\\minneigh{}",
                    lower_label     + " (\\%)",
                    mid_low_label   + " (\\%)",
                    mid_high_label  + " (\\%)",
                    total_label     + " (\\%)",
                    "All \\Isols{}"), 
                    file=file)
        print("\\\\", file=file)
        print("\\hline", file=file)
        print("\\hline", file=file)
        for (i,l,ml,mh,h) in zip(x,lower, mid_low, mid_high, total):
            clustered = sum((l,ml,mh))
            print("\\bf %d & %.1f & %.1f & %.1f & %.1f & %d" % (
                i,
                100 * l / h,
                100 * ml / h,
                100 * mh / h,
                100 * (h - clustered) / h,
                h),
            file=file)
            print("\\\\", file=file)
            print("\\hline", file=file)
        print("\\end{tabular}", file=file)
        print("\\caption{The percentage of all \\isols{} that ended up in clusters where their \\spec{} is the minority (%s), the dominant (%s), and make up the entirety of the cluster (%s) or were categorized as %s.}" % 
                (lower_label,mid_low_label,mid_high_label, total_label), file=file)
        print("\\label{tab:coverage:percent:total}", file=file)
        print("\\end{table}", file=file)
    # COVERAGE PERCENT OF CLUSTERED
    with open("coverage-percent-clustered.tex", "w") as file:
        print("\\begin{table}[ht!]", file=file)
        print("\\setlength{\\arrayrulewidth}{1pt}", file=file)
        print("\\newcolumntype{n}{>{\columncolor{cyan}}r}", file=file)
        print("\\newcolumntype{p}{>{\columncolor{cpgreen}}>{\color{white}}r}", file=file)
        print("\\newcolumntype{h}{>{\columncolor{cpgreengold}}r}", file=file)
        print("\\newcolumntype{m}{>{\columncolor{cpgold}}r}", file=file)
        print("\\centering", file=file)
        print("\\begin{tabular}{|r|m|h|p|r|}", file=file)
        print("\\hline", file=file)
        print("\\bf %s&\\bf %s&\\bf %s&\\bf %s&\\bf %s" %
                ("\\minneigh{}",
                    lower_label     + " (\\%)",
                    mid_low_label   + " (\\%)",
                    mid_high_label  + " (\\%)",
                    "Clustered \\Isols{}"),
                file=file)
        print("\\\\", file=file)
        print("\\hline", file=file)
        print("\\hline", file=file)
        for (i,l,ml,mh,h) in zip(x,lower, mid_low, mid_high, total):
            clustered = sum((l,ml,mh))
            print("\\bf %d&%.1f&%.1f&%.1f&%d" % (
                i,
                100 * l / clustered,
                100 * ml / clustered,
                100 * mh / clustered,
                clustered),
            file=file)
            print("\\\\", file=file)
            print("\\hline", file=file)
        print("\\end{tabular}", file=file)
        print("\\caption{The percent of \\isols{} clustered that ended up in a cluster where their \\spec{} is the minority (%s), the dominant (%s), and make up the entirety of the cluster (%s).}" % 
                (lower_label,mid_low_label,mid_high_label), file=file)
        print("\\label{tab:coverage:percent:clustered}", file=file)
        print("\\end{table}", file=file)
    # PERCENT CLUSTERED VS. UNCLUSTERED
    with open("coverage-percent-clustered-noise.tex", "w") as file:
        print("\\begin{table}[ht!]", file=file)
        print("\\setlength{\\arrayrulewidth}{1pt}", file=file)
        print("\\centering", file=file)
        print("\\begin{tabular}{|r|r|r|r|}", file=file)
        print("\\hline", file=file)
        print("\\bf %s&\\bf %s&\\bf %s&\\bf %s" %
                ("\\minneigh{}",
                    "Clustered" + " (\\%)",
                    "Noise"     + " (\\%)",
                    "All \\Isols{}"),
                file=file)
        print("\\\\", file=file)
        print("\\hline", file=file)
        print("\\hline", file=file)
        for (i,l,ml,mh,h) in zip(x,lower, mid_low, mid_high, total):
            clustered = sum((l,ml,mh))
            print("\\bf %d & %.1f & %.1f & %d" % (
                i,
                100 * clustered / h,
                100 * (h-clustered) / h,
                h),
            file=file)
            print("\\\\", file=file)
            print("\\hline", file=file)
        print("\\end{tabular}", file=file)
        print("\\caption{The percent of \\isols{} that \\dbscan{} placed in a cluster or determined to be noise.}", file=file)
        print("\\label{tab:coverage:percent:clustered:noise}", file=file)
        print("\\end{table}", file=file)
    # VALUES CLUSTERED VS. UNCLUSTERED
    with open("coverage-values-clustered-noise.tex", "w") as file:
        print("\\begin{table}[ht!]", file=file)
        print("\\setlength{\\arrayrulewidth}{1pt}", file=file)
        print("\\centering", file=file)
        print("\\begin{tabular}{|r|r|r|r|}", file=file)
        print("\\hline", file=file)
        print("\\bf %s&\\bf %s&\\bf %s&\\bf %s" %
                ("\\minneigh{}",
                    "Clustered",
                    "Noise",
                    "All \\Isols{}"),
                file=file)
        print("\\\\", file=file)
        print("\\hline", file=file)
        print("\\hline", file=file)
        for (i,l,ml,mh,h) in zip(x,lower, mid_low, mid_high, total):
            clustered = sum((l,ml,mh))
            print("\\bf %d & %d & %d & %d" % (
                i,
                clustered,
                h-clustered,
                h),
            file=file)
            print("\\\\", file=file)
            print("\\hline", file=file)
        print("\\end{tabular}", file=file)
        print("\\caption{The number of \\isols{} that \\dbscan{} placed in a cluster or determined to be noise.}", file=file)
        print("\\label{tab:coverage:values:clustered:noise}", file=file)
        print("\\end{table}", file=file)

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    stacked_lower    = [l + 0 for l in lower]
    stacked_mid_low  = [l + ml for l,ml in zip(lower, mid_low)]
    stacked_mid_high = [l + ml +  mh for l,ml,mh in zip(lower, mid_low, mid_high)]
    stacked_total    = [0 + t for t in total]

    lower_color = gold
    mid_low_color = green_gold
    mid_high_color = green
    total_color = 'magenta' # Hurts my eyes
    total_color = 'cyan'

    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    axes.fill_between(x, stacked_mid_high, stacked_total,    facecolor=total_color)
    axes.fill_between(x, stacked_mid_low,  stacked_mid_high, facecolor=mid_high_color)
    axes.fill_between(x, stacked_lower,    stacked_mid_low,  facecolor=mid_low_color)
    axes.fill_between(x, 0,                stacked_lower,    facecolor=lower_color)
    pyplot.scatter(x, stacked_lower,    label=lower_label,    c=lower_color, linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_mid_low,  label=mid_low_label,  c=mid_low_color, linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_mid_high, label=mid_high_label, c=mid_high_color, linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_total,    label=total_label,    c=total_color, linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_lower_upper_scatter_filled_stack(x, lower, upper, total, title,
        filename='cluster_scatter.png', x_low=None, x_high=None, y_low=None,
        x_label='Minimum Neighbors', y_label='Number of Isolates',
        lower_label='Minority', upper_label='Majority', total_label='Unclustered',
        legend_title='Isolate Clustering', flip_colors=True,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    print('lower:{}'.format(lower))
    print('upper:{}'.format(upper))
    print('total:{}'.format(total))

    assert(all(l + u <= t for l,u,t in zip(lower, upper, total)))

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    stacked_lower    = [l + 0 for l in lower]
    stacked_upper  = [l + u for l,u in zip(lower, upper)]
    stacked_total     = [0 + t for t in total]

    if flip_colors:
        lower_color = gold
        upper_color = green
    else:
        lower_color = green
        upper_color = gold
    total_color = 'magenta'
    total_color = 'cyan'

    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    axes.fill_between(x, 0,                stacked_total,   facecolor=total_color)
    axes.fill_between(x, 0,                stacked_lower, facecolor=lower_color)
    axes.fill_between(x, stacked_lower,    stacked_upper, facecolor=upper_color)
    pyplot.scatter(x, stacked_lower, label=lower_label, c=lower_color, linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_upper, label=upper_label, c=upper_color, linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_total, label=total_label, c=total_color, linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_cluster_scatter_filled_stack(x, clust, unclust, multiply, total, title,
        filename='cluster_scatter.png', x_low=None, x_high=None, y_low=None,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    # print('x       :{}'.format(x))
    # print('clust   :{}'.format(clust))
    # print('unclust :{}'.format(unclust))
    # print('multiply:{}'.format(multiply))
    # print('total   :{}'.format(total))

    assert(all(c + u == t for c,u,t in zip(clust, unclust, total)))

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    stacked_clust    = [c + 0 for c in clust]
    stacked_unclust  = [c + u for c,u in zip(clust, unclust)]
    stacked_multiply = [c - m for c,m in zip(clust, multiply)]
    stacked_total     = [0 + t for t in total]

    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel('Minimum Neighbors')
    pyplot.ylabel('Count')
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    #axes.fill_between(x, 0,                stacked_total,   facecolor='cyan')
    axes.fill_between(x, 0,                stacked_unclust, facecolor=green)
    axes.fill_between(x, stacked_clust,    stacked_unclust, facecolor=gold)
    #axes.fill_between(x, stacked_multiply, stacked_clust,   facecolor='magenta')
    pyplot.scatter(x, stacked_clust,    label='Clustered',   c=green,     linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_unclust,  label='Unclustered', c=gold,      linewidths=.5, s=marker_size)
    #pyplot.scatter(x, stacked_multiply, label='Multiply',    c='magenta', linewidths=.5, s=marker_size)
    pyplot.scatter(x, stacked_total,    label='Total',       c='cyan',   linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_cluster_scatter_filled_between(x, clust, unclust, multiply, total, title,
        filename='cluster_scatter.png', x_low=None, x_high=None, y_low=None,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    print('x       :{}'.format(x))
    print('clust   :{}'.format(clust))
    print('unclust :{}'.format(unclust))
    print('multiply:{}'.format(multiply))
    print('total   :{}'.format(total))

    assert(all(c + u == t for c,u,t in zip(clust, unclust, total)))

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel('Minimum Neighbors')
    pyplot.ylabel('Count')
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    axes.fill_between(x, 0, total,    facecolor='cyan')
    axes.fill_between(x, 0, unclust,  facecolor=gold)
    axes.fill_between(x, unclust, clust,    facecolor=green)
    axes.fill_between(x, 0, multiply, facecolor='magenta')
    pyplot.scatter(x, clust,    label='Clustered',   c=green,     linewidths=.5, s=marker_size)
    pyplot.scatter(x, unclust,  label='Unclustered', c=gold,      linewidths=.5, s=marker_size)
    #pyplot.scatter(x, multiply, label='Multiply',    c='magenta', linewidths=.5, s=marker_size)
    pyplot.scatter(x, total,    label='Total',       c='cyan',   linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_cluster_scatter_filled_below(x, clust, unclust, multiply, total, title,
        filename='cluster_scatter.png', x_low=None, x_high=None, y_low=None,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    print('x       :{}'.format(x))
    print('clust   :{}'.format(clust))
    print('unclust :{}'.format(unclust))
    print('multiply:{}'.format(multiply))
    print('total   :{}'.format(total))

    assert(all(c + u == t for c,u,t in zip(clust, unclust, total)))

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel('Minimum Neighbors')
    pyplot.ylabel('Count')
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    axes.fill_between(x, 0, total,    facecolor='cyan')
    axes.fill_between(x, 0, clust,    facecolor=green)
    axes.fill_between(x, 0, unclust,  facecolor=gold)
    axes.fill_between(x, 0, multiply, facecolor='magenta')
    pyplot.scatter(x, clust,    label='Clustered',   c=green,     linewidths=.5, s=marker_size)
    pyplot.scatter(x, unclust,  label='Unclustered', c=gold,      linewidths=.5, s=marker_size)
    #pyplot.scatter(x, multiply, label='Multiply',    c='magenta', linewidths=.5, s=marker_size)
    pyplot.scatter(x, total,    label='Total',       c='cyan',   linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_cluster_scatter_unfilled(x, clust, unclust, multiply, total, title,
        filename='cluster_scatter.png', x_low=None, x_high=None, y_low=None,
        y_high=None, marker_size=20):
    axes = pyplot.gca()
    print('x       :{}'.format(x))
    print('clust   :{}'.format(clust))
    print('unclust :{}'.format(unclust))
    print('multiply:{}'.format(multiply))
    print('total   :{}'.format(total))

    assert(all(c + u == t for c,u,t in zip(clust, unclust, total)))

    if x_low is None:
        x_low = min(x[0],x[-1]) - 1
    if x_high is None:
        x_high = max(x[0],x[-1]) + 1
    if x_low < x_high:
        axes.set_xlim(x_low, x_high)
    else:
        axes.set_xlim(x_high, x_low)
    y_delt = 500
    y_upper = max(total) + y_delt
    axes.set_ylim(0, y_upper)
    pyplot.xlabel('Minimum Neighbors')
    pyplot.ylabel('Count')
    pyplot.yticks(numpy.arange(0, y_upper, y_delt))
    pyplot.plot(x, clust,       c=green)
    pyplot.plot(x, unclust,   c=gold)
    #pyplot.plot(x, multiply,     c='magenta')
    pyplot.plot(x, total,           c='cyan')
    pyplot.scatter(x, clust,    label='Clustered',   c=green,     linewidths=.5, s=marker_size)
    pyplot.scatter(x, unclust,  label='Unclustered', c=gold,      linewidths=.5, s=marker_size)
    #pyplot.scatter(x, multiply, label='Multiply',    c='magenta', linewidths=.5, s=marker_size)
    pyplot.scatter(x, total,    label='Total',       c='cyan',   linewidths=.5, s=marker_size)
    pyplot.title(title)
    legend = pyplot.legend(loc=9, bbox_to_anchor=(0.5, -0.1), title='Isolate Clustering', ncol=4)
    pyplot.savefig(filename, bbox_extra_artists=(legend,), bbox_inches='tight');
    pyplot.close()

def plot_dist_clust(vals, title, filename='clust_dist.png'):
    pyplot.hist(vals, bins=[0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1], color=default_color)
    axes = pyplot.gca()
    axes.set_xlim(0, 1.05)
    pyplot.ylabel('Number of Clusters')
    pyplot.xlabel('Individual Cluster Purity')
    pyplot.xticks(numpy.arange(0, 1.04, .10))
    pyplot.title(title)
    pyplot.savefig(filename, bbox_inches='tight');
    pyplot.close()
    return
def plot_dist_size(vals, title, filename='size_dist.png'):
    binwidth = 10
    bins = numpy.arange(0, 380 + binwidth, binwidth)
    axes = pyplot.gca()
    axes.set_ylim(0, 80)
    pyplot.hist(vals, bins=bins, color=default_color)
    pyplot.title(title)
    #pyplot.xticks(bins)
    pyplot.ylabel('Number of Clusters')
    pyplot.xlabel('Size of Cluster')
    pyplot.savefig(filename, bbox_inches='tight');
    pyplot.close()
    return
def plot_dist_data(vals, title, filename='clust_dist_data.png'):
    pyplot.hist(vals, bins=[0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1], color=default_color)
    axes = pyplot.gca()
    axes.set_xlim(0, 1.05)
    axes.set_ylim(0, 1600)
    pyplot.ylabel('Number of Datapoints')
    pyplot.xlabel('Individual Cluster Purity')
    pyplot.xticks(numpy.arange(0, 1.04, .10))
    pyplot.title(title)
    pyplot.savefig(filename, bbox_inches='tight');
    pyplot.close()
    return

def main():
    print("Nothing implemented just yet.")
    return 0

if __name__ == '__main__':
    rtn = main()
    sys.exit(rtn)
