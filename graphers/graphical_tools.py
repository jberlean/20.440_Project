
from scipy.interpolate import spline

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import itertools
import jgraph as igraph





# This library consists of graphical tools that may be useful to visualize results

'''
Factory Methods
'''


def pairs(lst):
    return [(a,b) for i,a in enumerate(lst) for b in lst[i+1:]]


def find_vertices_edges(raw_info,set_vertices = False):
    if raw_info == []: return [],[]
    clones = [['A'+str(a) for a in c[0]] + ['B'+str(b) for b in c[1]] for c in raw_info]
    if set_vertices: vertices = set_vertices
    else: vertices = list(set([a for b in clones for a in b]))
    clones = [pairs(c) for c in clones]
    clones = [l for c in clones for l in c] # flatten list
    edges = [(vertices.index(c[0]),vertices.index(c[1])) for c in clones]
    return vertices,edges


def edge_frequencies(raw_info,freqs):
    if raw_info == []: return [],[]
    clones = [['A'+str(a) for a in c[0]] + ['B'+str(b) for b in c[1]] for c in raw_info]
    edge_freqs = [f for f,c in zip(freqs,clones) for _ in pairs(c)]
    edge_freqs = [np.log10(f) for f in edge_freqs] # apply function
    edge_freqs = [3*(f - min(edge_freqs)) + 1 for f in edge_freqs] # normalize
    return edge_freqs


def chain_frequencies(data):
    pass
#class MethodComparison:
#    def __init__(self, result_list, labels):
        

def graphical_frequencies(results,data,mod=True):
    actual_freq = []
    predicted_freq = []

    cells_per_well = data.metadata['cells_per_well_distribution_params']['cells_per_well']
    wells = data.metadata['num_wells']
    
    for c,f in zip(data.metadata['cells'],data.metadata['generated_data']['cell_frequencies']):
        if c in results['cells']:
            val = results['cell_frequencies'][results['cells'].index(c)]
            if mod: predicted_freq.append(1 - ((1-val)**(1./cells_per_well)))
            else: predicted_freq.append(val)
            actual_freq.append(f)
        else:
            pass#predicted_freq.append(0.)
    
    fig = plt.figure()
    ax = plt.gca()
    plt.scatter(predicted_freq,actual_freq)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('Predicted Frequency')
    plt.ylabel('Actual Frequency')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlim((0.0001,0.1))
    plt.ylim((0.0001,0.1))
    plt.show()

def graphical_data_summary(data,results=None):

    f, axarr = plt.subplots(1, 2)

    cells = data.metadata['cells']
    freqs = data.metadata['generated_data']['cell_frequencies']
    
    limits = (np.floor(np.log10(min(freqs))),np.ceil(np.log10(max(freqs))))
    x = np.logspace(limits[0], limits[1], 1*(limits[1]-limits[0])+1)
    hist = np.histogram(freqs, bins=x)

    for ind,i in enumerate(hist[0]):
        axarr[0].plot([x[ind],x[ind+1]],[i,i],'k--')
        print x[ind],x[ind+1]

    coverage = []
    for ind,r in enumerate(results):
        correct_freqs = [freqs[cells.index(p)] for p in r['cells'] if p in cells]
        coverage.append(np.histogram(correct_freqs, bins=x))
        #print 'Coverage:',coverage[-1]
        #print 'X:',x
        axarr[0].bar(x[:-1] + (ind+1)*x[:-1]*(10**0.25),coverage[-1][0],x[:-1]*(10**0.25))
    
            
    

    axarr[0].set_xscale('log')
    axarr[0].set_xlim([x[0],x[-1]])
    
    plt.show(block=False)
    raw_input('Press enter to close...')

    #axarr[0, 0].plot(x, y)
    
    
def graphical_network(data,results1 = {'cells':[]},results2 = {'cells':[]},hide_labels = True):
    
    vertices,edges = find_vertices_edges(data.metadata['cells'])
    edge_width = edge_frequencies(data.metadata['cells'],data.metadata['generated_data']['cell_frequencies'])
    vertice_size = chain_frequencies(data)
    _,edges_1 = find_vertices_edges(results1['cells'],set_vertices = vertices)
    _,edges_2 = find_vertices_edges(results2['cells'],set_vertices = vertices)
    #edges.sort(key=lambda x: x[0])
    #edges_1.sort(key=lambda x: x[0])
    #edges_2.sort(key=lambda x: x[0])

    # generate edge colors corresponding to prediction success
    edge_colors = []
    for edge in edges:
        set1,set2 = edge in edges_1,edge in edges_2
        if all((set1,set2)): edge_colors.append('black')
        elif set1: edge_colors.append('green')
        elif set2: edge_colors.append('yellow')
        else: edge_colors.append("grey")          
            
    visual_style = {}
    visual_style["vertex_size"] = 20
    visual_style["vertex_color"] = ['red' if v[0] == 'A' else 'blue' for v in vertices]
    visual_style["edge_color"] = edge_colors
    visual_style["edge_width"] = 15
    visual_style["vertex_size"] = 7
    visual_style["edge_width"] = edge_width
    #visual_style["vertex_label"] = g.vs["name"]
    #visual_style["edge_width"] = [1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]
    #visual_style["layout"] = layout
    #visual_style["bbox"] = (300, 300)
    #visual_style["margin"] = 20
    
    if hide_labels: vertices = ['' for v in vertices]

    g = igraph.Graph(vertex_attrs={"label": vertices}, edges=edges, directed=False)

    igraph.plot(g,**visual_style)

    
def graphical_auroc(results1,results2,seqdata):

    def compute_roc_curve(guess, threshold, max_tpr, max_fpr):
        data_sorted = sorted(zip(guess, threshold), key=lambda v: v[1], reverse=True)

        x,y = [0.],[0.]
        last_t = max(threshold) + 1
        for g,t in data_sorted:
            if g in hits:
                dx, dy = 0, 1./max_tpr
            else:
                dx, dy = 1./max_fpr, 0
            if last_t != t:
                x.append(x[-1]+dx)
                y.append(y[-1]+dy)
            else:
                x[-1] += dx
                y[-1] += dy
            last_t = t
        x.append(1.0)
        y.append(1.0)

        return x, y
 
    obs_alphas, obs_betas = map(lambda l: set(sum(l, [])), zip(*seqdata.well_data))

    hits = set(filter(lambda c: all([a in obs_alphas for a in c[0]] + [b in obs_betas for b in c[1]]), seqdata.metadata['cells']) )
    guess1,guess2 = results1['cells'],results2['cells']
    thresh1,thresh2 = results1['threshold'],results2['threshold']

    pos_count = len(hits)
    #neg_count = max(len(set(guess1)-set(hits)), len(set(guess2)-set(hits)))
    neg_count = len(obs_alphas)*len(obs_betas)*(1 + len(obs_alphas) + len(obs_betas)) - len(hits)

    #pos_count1,pos_count2 = len([g for g in guess1 if g in hits]),len([g for g in guess2 if g in hits])
    #neg_count1,neg_count2 = len(guess1)-pos_count1,len(guess2)-pos_count2
    #pos_count = max((pos_count1,pos_count2))
    
    #lims = np.arange(1.,0.,-0.001)
    
    #x1,y1 = [0],[0]
    #x2,y2 = [0],[0]
    
    #for lim in lims:
#        # check for set 1
#        tpr,fpr = 0.,0.
#        for g,t in zip(guess1,thresh1):
#            if t >= lim:
#                if g in hits: tpr += 1
#                else: fpr += 1
#        x1.append(fpr/neg_count1)
#        y1.append(tpr/pos_count)
#        
#        # check for set 2
#        tpr,fpr = 0.,0.
#        for g,t in zip(guess2,thresh2):
#            if t >= lim:
#                if g in hits: tpr += 1
#                else: fpr += 1
#        x2.append(fpr/neg_count2)
#        y2.append(tpr/pos_count)

#    x1,x2 = [x/max(x1) for x in x1],[x/max(x2) for x in x2]

    # Compute data for set 1
    x1,y1=compute_roc_curve(guess1, thresh1, pos_count, neg_count)
    x2,y2=compute_roc_curve(guess2, thresh2, pos_count, neg_count)

    print min(filter(lambda v: v!=0, thresh1))

    
    
#    print 'Start 1:'
#    for x,y in zip(x1,y1)[-10:]:
#        print x,y
#    print 'Start 2:'
#    for x,y in zip(x2,y2)[-10:]:
#        print x,y
    
    
    print 'Method 1 auROC: {}'.format(np.trapz(y1,x1))    
    print 'Method 2 auROC: {}'.format(np.trapz(y2,x2))
    
    plt.ion()

    plt.figure(1)

    x1_count = [v*neg_count for v in x1[:-1]]
    y1_count = [v*pos_count for v in y1[:-1]]
    x2_count = [v*neg_count for v in x2[:-1]]
    y2_count = [v*pos_count for v in y2[:-1]]
    plt.plot(x1_count,y1_count)
    plt.plot(x2_count,y2_count)
    plt.plot([0,1./9*max(y1_count+y2_count)], [0, max(y1_count+y2_count)], 'k--')
    plt.plot([0,1./2*max(y1_count+y2_count)], [0, max(y1_count+y2_count)], 'k--')
    plt.plot([0,max(y1_count+y2_count)], [0, max(y1_count+y2_count)], 'k--')
    plt.xlim((0, max(x1_count+x2_count)))
    plt.ylim((0, max(y1_count+y2_count)))
    plt.xlabel('False positive count')
    plt.ylabel('True positive count')
    plt.legend(loc='best')

    plt.figure(2)
   
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.xlim((0.,1.))
    plt.ylim((0.,1.))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')


    plt.show()
                    
                
    
    







