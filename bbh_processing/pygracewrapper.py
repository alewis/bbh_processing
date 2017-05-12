from PyGrace.grace import Grace
from PyGrace.colors import ColorBrewerScheme
from subprocess import call

def set_graphtype(graph, graphtype):
    if graphtype=="linxy":
        graph.linxy()
    elif graphtype=="linx":
        graph.linx()
    elif graphtype=="liny":
        graph.liny()
    elif graphtype=="logx":
        graph.logx()
    elif graphtype=="logy":
        graph.logy()
    elif graphtype=="logxy":
        graph.logxy()
    elif graphtype=="logy":
        graph.logy()
    else: 
        raise ValueError("Invalid graphtype: "+graphtype)

def from_axis(ax, data, labels=None, 
              graphtype="linxy",
              dirname="grace_sandbox/", pltname="testplot.agr",
              symbolsize=0.,
              display=True):
    grace = Grace(colors=ColorBrewerScheme('Paired'))
    graph = grace.add_graph()
    
    title = str(ax.get_title())
    xlabel = str(ax.get_xlabel())
    ylabel = str(ax.get_ylabel())
    plot_array(data, labels=labels, title=title, xlabel=xlabel, ylabel=ylabel, 
               graphtype=graphtype, pltname=pltname, 
              symbolsize=symbolsize, display=display)

def display(dirname="grace_sandbox/", pltname="testplot.agr"):
    outname = dirname+pltname
    call(["xmgrace", outname]) 

def plot_array(data, 
               labels=None,
               title="Title", xlabel="Array 1", ylabel="Array 2",
               graphtype="linxy",
               dirname="grace_sandbox/", pltname="testplot.agr",
               symbolsize=0.,
               display=True):
    
    grace = Grace(colors=ColorBrewerScheme('Paired'))
    graph = grace.add_graph()
    graph.title.text = title
    graph.xaxis.label.text = xlabel
    graph.yaxis.label.text = ylabel

    if labels is not None:
        dsets = [graph.add_dataset(d, legend=l) for d, l in zip(data, labels)]
    else:
        dsets = [graph.add_dataset(d) for d in data]
    
    graph.set_different_colors()
    graph.set_different_symbols()
    for dset in dsets:
        dset.symbol.size=symbolsize
    
    #graph.set_world_to_limits()
    set_graphtype(graph, graphtype)
    graph.autoscale()

    outname = dirname+pltname
    grace.write_file(outname)
    if display:
       call(["xmgrace", outname]) 


