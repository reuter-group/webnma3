# DESCRIPTION: Common functions for ploting and data file output
# CONTRIBUTORS: Dandan Xue (dandan.xue@uib.no)

from os.path import join
from collections.abc import Iterable
import matplotlib.pyplot as plt

DIG = 3

def plot_and_record(xtuple, ytuple, title='', tar_dir='',
                    name='fig.png', style='default', **args):
    xs, xlabel = xtuple
    ys, ylabel = ytuple

    # use 1-indexed x-axis
    xs = [x+1 for x in xs]
    
    if style == 'bar':
        plt.bar(xs, ys)
        plt.xticks([e for e in xs if e%5==0])

    elif style == 'heatmap':    
        ys = [y+1 for y in ys] # use 1-indexed y-axis
        plt.contourf(xs, ys,
                     args['z'],
                     10,  
                     cmap = plt.cm.bwr,
                     vmax = 1,
                     vmin = -1,
        )
        plt.colorbar()        
    elif style == 'fluc_ss':
        plt.plot(xs, ys, lw=0.7)
        legend = {'E':'β-sheet', 'G':'Helix'}
        color = {'E': 'darkturquoise', 'G': 'red'}
        
        for ss,spans in args['ss_dict'].items():
            leg_ss  = legend[ss]    
            for span in spans:
                if span[1]-span[0] >= 3:  # ignore very small spans
                    plt.axvspan(span[0], span[1], facecolor=color[ss],\
                                alpha=0.5, label=leg_ss)
                    leg_ss = None
        if 'pdbs' in args:
            # TODO fix
            plt.legend(tuple(args['pdbs'] + ['Helix','beta']))
        else:    
            plt.legend()
    else:
        plt.plot(xs, ys, lw=0.7)
        
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    fig_name = join(tar_dir, name)
    plt.savefig(fig_name, format='png', dpi=200)
    plt.gcf().clear()

    if 'data_filename' in args:
        data_filename = args.pop('data_filename')
    else:
        data_filename = name[:-4] + '.txt'
    record(xs, ys, tar_dir, data_filename, style, **args)

    
def record(xs, ys, tar_dir, data_filename, style='default', **args):
    data_filename = join(tar_dir, data_filename)
    with open(data_filename,'w') as f:
        if style != 'heatmap':
            if 'header' in args:
                f.write(args['header']+'\n') 
            for x,y in zip(xs,ys):
                if isinstance(y, Iterable):
                    y_dig = [round(n, DIG) for n in y]
                    y_str = '\t'.join(map(str, y_dig)) 
                    f.write('{}\t{}\n'.format(x,y_str))
                else:                
                    f.write('{}\t{}\n'.format(x,round(y,DIG)))
        else:
            one_line = lambda l: ' '.join([str(round(e,DIG)) for e in l])
            f.writelines([one_line(l) + '\n' for l in args['z']])
    
