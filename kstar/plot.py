from encodings.idna import dots
import string
import numpy as np
import fpdf
import pandas as pd
import os
import fpdf

from kstar import helpers

#plotting packages
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns





class OrientationError(Exception):
    def __init__(self, message = "Orientation Invalid. Valid Orientations are : ", valid_orientations = ['left', 'right', 'top', 'bottom']):
        self.message = message + ', '.join(valid_orientations)
    def __str__(self):
        return self.message
    

def plot_jaci_between_samples(evidence,samples,title='', ax = None, annot = True, cluster = False, **kwargs):
    """
    This function creates a heatmap based on jaccard index of phosphopeptide identities

    Parameters
    ----------
    evidence : pandas DataFrame
        binarized dataframe indicating presence/absence of phosphopeptides in samples
    samples : list
        list of sample column names in the evidence dataframe to use for calculating jaccard index
    title : str, optional
        title of the heatmap
    ax : matplotlib Axes instance, optional
        axes to plot heatmap on. If None, new figure and axes created
    annot : bool, optional
        whether to annotate heatmap with jaccard index values
    cluster : bool, optional
        whether to cluster samples based on jaccard index before plotting heatmap
    **kwargs : additional keyword arguments
        additional keyword arguments to pass to seaborn heatmap function
    """
    jaccard_matrix = helpers.jaci_matrix_between_samples(evidence,samples)

    if cluster:
        #cluster row and columns of jaccard matrix
        linkage_matrix = linkage(jaccard_matrix, method='average', metric='euclidean')
        dendro = dendrogram(linkage_matrix, no_plot=True)
        clustered_order = dendro['leaves']
        jaccard_matrix = jaccard_matrix.iloc[clustered_order, clustered_order]

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))

    sns.heatmap(jaccard_matrix, annot=annot,fmt='.2f', cbar_kws={'label': 'Jaccard Index'}, vmin=0, vmax=1, ax=ax, **kwargs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 35, ha = 'right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
    ax.set_title(title)  

    #adjust colorbar font size
    cbar = ax.collections[0].colorbar  # Get the colorbar
    cbar.ax.tick_params(labelsize=9) 
    return ax

class DotPlot:
    """
    The DotPlot class is used for plotting dotplots, with the option to add clustering and context plots.
    The size of the dots based on the values dataframe, where the size of the dot is the area of the value * dotsize
           
    Parameters
    ----------
    values: pandas DataFrame instance
        values to plot 
    fpr : pandas DataFrame instance 
        false positive rates associated with values being plotted
    alpha: float, optional
        fpr value that defines the significance cutoff to use when plt
        default : 0.05
    inclusive_alpha: boolean
        whether to include the alpha (significance <= alpha), or not (significance < alpha).
        default: True
    binary_sig: boolean, optional
        indicates whether to plot fpr with binary significance or as a change color hue
        default : True
    dotsize : float, optional
        multiplier to use for scaling size of dots
    colormap : dict, optional
        maps color values to actual color to use in plotting
        default : {0: '#6b838f', 1: '#FF3300'}
    labelmap = 
        maps labels of colors, default is to indicate FPR cutoff in legend
        default : None
    facecolor : color, optional
        Background color of dotplot
        default : 'white'
    legend_title : str, optional
        Legend Title for dot sizes, default is `p-value'
    size_number : int, optional 
        Number of dots to attempt to generate for dot size legend
    size_color : color, optional
        Size Legend Color to use 
    color_title : str, optional
        Legend Title for the Color Legend
    markersize : int, optional
        Size of dots for Color Legend
    legend_distance : int, optional
        relative distance to place legends 
    figsize : tuple, optional 
        size of dotplot figure
    title : str, optional
        Title of dotplot
    xlabel : bool, optional
        Show xlabel on graph if True
    ylabel : bool, optional
        Show ylabel on graph if True
    x_label_dict: dict, optional
        Mapping dictionary of labels as they appear in values dataframe (keys) to how they should appear on plot (values)
    kinase_dict: dict, optional
        Mapping dictionary of kinase names as they appear in values dataframe (keys) to how they should appear on plot (values)
    
    Attributes
    ----------
    values: pandas dataframe
        a copy of the original values dataframe
    fpr: pandas dataframe
        a copy of the original fpr dataframe
    alpha: float
        cutoff used for significance, default 0.05
    inclusive_alpha: boolean
        whether to include the alpha (significance <= alpha), or not (significance < alpha)
    significance: pandas dataframe
        indicates whether a particular kinases activity is significant, where fpr <= alpha is significant, otherwise it is insignificant
    colors: pandas dataframe
        dataframe indicating the color to use when plotting: either a copy of the fpr or significance dataframe
    binary_sig: boolean
        indicates whether coloring will be done based on binary significance or fpr values. Default True
    labelmap: dict
        indicates how to label each significance color
    figsize: tuple
        size of the outputted figure, which is overridden if axes is provided for dotplot
    title: string
        title of the dotplot
    xlabel: boolean
        indicates whether to plot x-axis labels
    ylabel: boolean
        indicates whether to plot y-axis labels
    colormap: dict
        colors to be used when plotting
    facecolor: string
        background color of dotplot
    """
    
    
    def __init__(self, values, fpr, alpha = 0.05, inclusive_alpha = True,
                 binary_sig = True, dotsize = 5, 
                 colormap={0: '#6b838f', 1: '#FF3300'}, facecolor = 'white',
                 legend_title = '-log10(p-value)', size_number = 5, size_color = 'gray', 
                 color_title = 'Significant', markersize = 10, 
                 legend_distance = 1.0, figsize = (4,8), title = None,
                 xlabel = True, ylabel = True, x_label_dict = None, kinase_dict = None):

        #detect if values are already log-transformed. If not, transform them
        if np.all((values >= 0) & (values <= 1)):
            values = -np.log10(values)

        self.values = values.copy()
        self.fpr = fpr.copy()
        #make sure that fpr dataframe has the same index as values dataframe. If not, reindex
        self.fpr = self.fpr.loc[self.values.index,self.values.columns]
        self.alpha = alpha
        #create binary dataframe that indicates significance based on provided fpr cutoff.
        #if inclusive_alpha:
        #    self.significance = (self.fpr <= alpha) * 1
        #else:
        #    self.significance = (self.fpr < alpha) * 1
        #Assign either fpr or significance to colors dataframe based on 
        self.binary_sig = binary_sig
        self.inclusive_alpha = inclusive_alpha
        self.set_colors()
      
        
        self.figsize =  figsize
        self.title = title
        self.xlabel = xlabel
        self.ylabel= ylabel

        self.colormap = colormap
        
        self.facecolor = facecolor

        self.dotsize = dotsize
        
        self.legend_title = legend_title
        self.size_number = size_number
        self.size_color = size_color
        self.markersize = markersize

        self.color_title = color_title
        self.legend_distance = legend_distance

        self.multiplier = 10
        self.offset = 5

        self.columns = self.set_column_labels(x_label_dict)
        self.index = self.set_index_labels(kinase_dict)
    
    def set_values(self, values):
        self.values = values


    def set_colors(self, labelmap = None):
        """
        Set colors for the plot based on significance or false positive rate.
        """
        #create binary dataframe that indicates significance based on provided fpr cutoff.
        if self.inclusive_alpha:
            self.significance = (self.fpr <= self.alpha) * 1
        else:
            self.significance = (self.fpr < self.alpha) * 1

        #define the colors to use when plotting (either gradient or binary)
        if self.binary_sig:
            #set colors to significance dataframe
            self.colors = self.significance
            if labelmap is None:
                if self.inclusive_alpha:
                    self.labelmap = {0: 'FPR > %0.2f'%(self.alpha), 1:'FPR <= %0.2f'%(self.alpha)}
                else: 
                    self.labelmap = {0: 'FPR >= %0.2f'%(self.alpha), 1:'FPR < %0.2f'%(self.alpha)}
            else:
                self.labelmap = labelmap

        else:
            self.colors = self.fpr

    def set_column_labels(self, x_label_dict):
        self.column_labels = list(self.values.columns)

        if x_label_dict is None: #just strip the data: string
            self.x_label_dict = {}
            
            #build an x_label_dict 
            for col in self.column_labels:
                self.x_label_dict[col] = col.replace('data:','')
            self.column_labels = [x.replace('data:','') for x in self.column_labels]

        else:
            #check that the label dictionary keys matches the columns
            labels = x_label_dict.keys()
            if set(labels) != set(self.column_labels):
                raise ValueError("The x_label_dict must have the same elements as the value columns")
            else:
                label_arr = []
                for col in self.column_labels:
                    label_arr.append(x_label_dict[col])
            self.column_labels = label_arr
            self.x_label_dict = x_label_dict
            
    def set_index_labels(self, kinase_dict):
        self.index_labels = list(self.values.index)
        if kinase_dict is None:
            self.kinase_dict = kinase_dict
        elif isinstance(kinase_dict, dict):
            #if custom dictionary is provided, make sure the appropriate elements are found inside it (needs to be at least
            names = kinase_dict.keys()
            if not set(self.index_labels).issubset(set(names)):
                raise ValueError("The kinase_dict must contain at least all the kinases found in values")
            else:
                label_arr = []
                for index in self.index_labels:
                    label_arr.append(kinase_dict[index])
            
                self.index_labels = label_arr
                self.kinase_dict = kinase_dict
        else:
            raise TypeError("If wanting to do a custom naming system, a custom dictionary must be provided in the 'kinase_dict' parameter")
        

    def dotplot(self, ax = None, orientation = 'left', size_legend = True, color_legend = True, max_size = None, **kwargs):
        """
        Generates the dotplot plot, where size is determined by values dataframe and color is determined by significant dataframe
        
        Parameters
        -----------
        ax : matplotlib Axes instance, optional
            axes dotplot will be plotted on. If None then new plot generated
        orientation : str, optional
            orientation to place legends, either 'left' or 'right'
        size_legend : bool, optional
            whether to include size legend (indicates meaning of dot size/activity)
        color_legend : bool, optional
            whether to include color legend (indicates significance)
        max_size : int, optional
            maximum size value to use when generating size legend. If None, automatic legend generated
        
        Returns
        -------
        ax : matplotlib Axes instance
            Axes containing the dotplot
        """
        valid_orientations = ['left', 'right']
        if orientation not in valid_orientations:
            raise OrientationError(valid_orientations = valid_orientations)
            
        if not ax:
            fig, ax = plt.subplots(figsize=self.figsize)
        ax.set_facecolor(self.facecolor)
        ax.set_title(self.title)
        
        # Transform Data
        columns = list(self.values.columns)
        self.values['row_index'] = np.arange(len(self.values)) * self.multiplier + self.offset
        self.colors['row_index'] = np.arange(len(self.colors)) * self.multiplier + self.offset
    
        melt = self.values.melt(id_vars = 'row_index')
        self.values.drop(columns = ['row_index'], inplace = True)
        melt['var'] = melt.apply(lambda row : columns.index(row.iloc[1]) * self.multiplier + self.offset, axis = 1)
        
        melt_color = self.colors.melt(id_vars = 'row_index')
        melt_color['var'] = melt_color.apply(lambda row : columns.index(row.iloc[1]) * self.multiplier + self.offset, axis = 1)
        self.colors.drop(columns = ['row_index'], inplace = True)

        # Plot Data
        x = melt['var']
        y = melt['row_index'][::-1]    #needs to be done in reverse order to maintain order in the dataframe
        
        
        s = melt.value * self.dotsize
        
        #check to see if more than 2 values are given (fprs). Otherwise get color based on binary significance
        if self.binary_sig:
            #get color for each datapoint based on significance
            melt_color['color'] = [self.colormap.get(l,'black') for l in melt_color.value]
        else:
            cmap = LinearSegmentedColormap.from_list("sig_cmap", [self.colormap[0], self.colormap[1]])
            norm = Normalize(vmin=0, vmax=2, clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
            #replace 0 with 0.01 to avoid log10 errors, transform the fprs with a log transform
            melt_color.replace(0, 0.01, inplace=True)
            melt_color.value = -np.log10(melt_color.value)
            #get color for each datapoint based on fpr value
            melt_color['color'] = [mapper.to_rgba(l) for l in melt_color.value]

            
        c = melt_color['color']
        scatter = ax.scatter(x, y, c=c, s=s, **kwargs)
        
        # Add Color Legend
        if color_legend:
            if self.binary_sig:
                #create the legend
                color_legend = []
                for color_key in self.colormap.keys():
                    color_legend.append(
                        Line2D([0], [0], marker='o', color='w', label=self.labelmap[color_key],
                                markerfacecolor= self.colormap[color_key], markersize=self.markersize),
                    )     
                legend1 = ax.legend(handles=color_legend, loc=f'upper {orientation}', bbox_to_anchor=(self.legend_distance,1), title = self.color_title)  
            else:
                #choose which values to show in the legend
                legend_vals = [1, 0.5, 0.05, 0.01]
                legend_color = [mapper.to_rgba(-np.log10(val)) for val in legend_vals]
                #create the legend 
                color_legend = []
                for i in range(len(legend_vals)):
                    color_legend.append(Line2D([0], [0], marker='o', color='w', label=str(legend_vals[i]),
                                markerfacecolor= legend_color[i], markersize=self.markersize))
                legend1 = ax.legend(handles=color_legend, loc=f'upper {orientation}', bbox_to_anchor=(self.legend_distance,1), title = 'FPR')  
            legend1.set_clip_on(False)
            ax.add_artist(legend1)
            


        # Add Size Legend
        if size_legend:
            #check to see if max pval parameter was given: if so, use to create custom legend
            if max_size is not None:
                s_label = np.arange(max_size/self.size_number,max_size+1,max_size/self.size_number).astype(int)
                dsize = [s*self.dotsize for s in s_label]
                legend_elements = []
                for element, s in zip(s_label, dsize):
                    legend_elements.append(Line2D([0],[0], marker='o', color = 'w', markersize = s**0.5, markerfacecolor = self.size_color, label = element))
                legend2 = ax.legend(handles = legend_elements, loc = f'lower {orientation}', title = self.legend_title, bbox_to_anchor=(self.legend_distance,0))        
            else:
                kw = dict(prop="sizes", num=self.size_number, color=self.size_color, func=lambda s: s/self.dotsize) 
                legend2 = ax.legend(*scatter.legend_elements(**kw),
                        loc=f'lower {orientation}', title=self.legend_title, bbox_to_anchor=(self.legend_distance,0)) 
            #ax.add_artist(legend2)

        
        # Add Additional Plotting Information
        ax.tick_params(axis = 'x', rotation = 90)
        ax.yaxis.set_ticks(np.arange(len(self.values)) * self.multiplier + self.offset)
        ax.xaxis.set_ticks(np.arange(len(columns)) * self.multiplier + self.offset)
        
        # set column labels in case values has changed
        self.set_column_labels(self.x_label_dict)
        ax.set_xticklabels(self.column_labels)
        ax.set_yticklabels(self.index_labels[::-1])
        #adjust x and y scale so that data is always equally spaced
        ax.set_ylim([0,len(self.values)*self.multiplier])
        ax.set_xlim([0,len(columns)*self.multiplier])
        
        if not self.xlabel:
            ax.axes.xaxis.set_visible(False)
        if not self.ylabel:
            ax.axes.yaxis.set_visible(False)
        return ax 
    
    def cluster(self, ax, method='single', metric='euclidean', orientation = 'top', color_threshold = -np.inf):
        """
        Performs hierarchical clustering on data and plots result to provided Axes. 
        result and significant dataframes are ordered according to clustering
        
        Parameters
        ----------
        ax : matplotlib Axes instance
            Axes to plot dendogram to
        
        method : str, optional
            The linkage algorithm to use.
        metric : str or function, optional
            The distance metric to use in the case that y is a collection of observation vectors; 
            ignored otherwise. See the pdist function for a list of valid distance metrics. A custom distance function can also be used.
        
        orientation : str, optional
            The direction to plot the dendrogram, which can be any of the following strings:
            'top': Plots the root at the top, and plot descendent links going downwards. (default).
            'bottom': Plots the root at the bottom, and plot descendent links going upwards.
            'left': Plots the root at the left, and plot descendent links going right.
            'right': Plots the root at the right, and plot descendent links going left.
        """
        if orientation in ['left', 'right']:
            row_linkage = linkage(self.values, method = method, metric = metric)
            den_row = dendrogram(row_linkage, 
                        ax = ax, 
                        orientation = orientation, 
                        labels = list(self.values.index), 
                        color_threshold = color_threshold, 
                        above_threshold_color = 'black', 
                        no_labels = True, 
                        show_leaf_counts = False) 
            self.values = self.values.iloc[den_row['leaves']].copy()
            self.colors = self.colors.iloc[den_row['leaves']].copy()
            self.set_index_labels(self.kinase_dict)

        
        elif orientation in ['top', 'bottom']:
            col_linkage = linkage(self.values.T, method=method, metric = metric)
            den_col = dendrogram(col_linkage, 
                            ax = ax, 
                            orientation = orientation, 
                            labels = list(self.values.columns), 
                            color_threshold = color_threshold, 
                            above_threshold_color = 'black', 
                            no_labels = True, 
                            show_leaf_counts = False)
            self.values = self.values.iloc[:, den_col['leaves']].copy()
            self.colors = self.colors.iloc[:,den_col['leaves']].copy()
            self.set_column_labels(self.x_label_dict)

        else:
            raise OrientationError()
            
        ax.tick_params(axis='both', which='both', length=0)
        
    def drop_kinases_with_no_significance(self):
        """
        Drop kinases from the values dataframe (inplace) when plotting if they are never observed as significant
        
        """

        kinase_list = self.significance[self.significance.sum(axis=1) ==0].index.values
        #check to make sure kinase_list only contains kinases currently in values dataframe
        kinase_list = [kin for kin in kinase_list if kin in self.index_labels]
        #remove kinases
        self.drop_kinases(kinase_list)


    def drop_kinases(self, kinase_list):
        """
        Given a list of kinases, drop these from the dot.values dataframe in all future plotting of this object. Removal 
        is in place

        Parameters
        ----------
        kinase_list: list
            list of kinase names to remove

        """ 
        #check to make sure kinase_list only contains kinases currently in values dataframe
        kinase_list = [kin for kin in kinase_list if kin in self.index_labels]
        
        self.values.drop(index=kinase_list, inplace=True)
        self.colors.drop(index = kinase_list, inplace=True)

        #update index_labels property as well
        for kin in kinase_list:
            if self.kinase_dict is None:
                self.index_labels.remove(kin)
            else:
                self.index_labels.remove(self.kinase_dict[kin])
        
        
    def context(self, ax, info, id_column, context_columns, dotsize = 200, markersize = 20, orientation = 'left', color_palette='colorblind', margin = 0.2, make_legend = True, **kwargs):
        """
        Context plot is generated and returned. The context plot contains the categorical data used for describing the data.
        
        Parameters
        ----------
        ax : maptlotlib axis
            where to map subtype information to
        info : pandas df
            Dataframe where context information is pulled from
        id_column: str
            Column used to map the subtype information to
        context_columns : list
            list of columns to pull context informaiton from
        dotsize : int, optional
            size of context dots
        markersize: int, optional
            size of legend markers
        orientation : str, optional
            orientation to plot context plots to - determines where legends are placed
            options : left, right, top, bottom
        color_palette : str, optional
            seaborn color palette to use  
        margin: float, optional
            margin  
        make_legend : bool, optional
            whether to create legend for context colors
        """
        
        orientation_values = {
            'left' : -1,
            'right' : 1,
            'bottom' : -1,
            'top' : 1
        }
        
        if orientation in [ 'left', 'right']:
            index = list(self.values.index)[::-1]  #reverse order for left/right orientation
        elif orientation in ['top', 'bottom']:
            index = list(self.values.columns)
        else: 
            raise OrientationError
        
        #record the number of different context types to include
        num_context = len(context_columns)
        melted = info[[id_column] + context_columns].melt(id_vars=id_column)

        #weird issue with melt function here, where for one datset it provides the context column names in 0 column rather than 'variable'. Rename for now.
        if 0 in melted.columns:
            melted.rename(columns = {0: 'variable'}, inplace = True)
        melted['var'] = melted.apply(lambda row : index.index(row.iloc[0]) * self.multiplier + self.offset, axis = 1)
        color_labels = melted['value'].unique()
        rgb_values = sns.color_palette(color_palette, len(color_labels))
        color_map = dict(zip(color_labels, rgb_values))
        
        if orientation in ['left', 'right']:
            ax.scatter(x = melted['variable'], y = melted['var'],c = melted['value'].map(color_map), s = dotsize, **kwargs)
            ax.tick_params(axis = 'x', rotation = 90)
            ax.axes.get_yaxis().set_visible(False)
            ax.margins(margin, 0.05)
        elif orientation in ['top', 'bottom']:
            ax.scatter(x = melted['var'], y = melted['variable'], c = melted['value'].map(color_map), s = dotsize, **kwargs)
            ax.axes.get_xaxis().set_visible(False)
            ax.margins(0.05, margin)
            ax.set_ylim([-0.5,num_context+0.5-1])
        
           
        total = len(melted['value'].unique()) + len(info.columns)-1
        running_total = 0
        
        # Add legends
        if make_legend:
            for col in context_columns:
                ids = info[col].unique()
                sig_legend = []
                for label in ids:
                    color = color_map[label]
                    sig_legend.append(
                        Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor=color,markersize=markersize))
                    if orientation in ['left', 'right']:
                        leg = ax.legend(
                            handles=sig_legend, 
                            bbox_to_anchor=(orientation_values[orientation],1-running_total/total), 
                            title=col)
                    elif orientation in ['top', 'bottom']:
                        leg = ax.legend(
                            handles=sig_legend, 
                            bbox_to_anchor=(running_total/total, orientation_values[orientation]),
                            loc='lower left',
                            title=col)
                leg.set_clip_on(False)
                ax.add_artist(leg)
    
                running_total += len(ids) + 1
                
    def evidence_count(self, ax, binary_evidence, plot_type = 'bars', phospho_type = None, dot_size = 1, include_recommendations = False,
                      ideal_min = None, recommended_min = None, dot_colors = None, bar_line_colors = None):
        """
        Add bars to dotplot indicating the total number of sites used as evidence in activity calculation

        Parameters
        ----------
        ax: axes object
            where to plot the bars
        binary_evidence: pandas dataframe
            binarized dataframe produced during activity calculation (threshold applied to original experiment)
        """
        palette = sns.color_palette('colorblind')
        
        #make sure binary evidence contains unique sites
        binary_evidence = binary_evidence.drop_duplicates()

        #get the number of sites used as evidence
        total_num_sites = binary_evidence.shape[0]
        num_sites_in_sample = binary_evidence[self.x_label_dict.keys()].sum()

        #get correct order of samples
        order = self.values.columns
        xticks = [tick * self.multiplier + self.offset for tick in range(len(order))]
        num_sites_in_sample = num_sites_in_sample[order]
        
        #add recommended labels
        if include_recommendations:
            #if minimums not indicated, use our recommendations
            if phospho_type == 'Y':
                if ideal_min is None:
                    ideal_min = 50
                if recommended_min is None:
                    recommended_min = 25
            elif phospho_type == 'ST':
                if ideal_min is None:
                    ideal_min = 1000
                if recommended_min is None:
                    recommended_min = 250
            else:
                raise ValueError('Please indicate which phospho_type this plot is for to get appropriate recommendations about the recommended evidence sizes (phospho_type = "Y")')

            #calculate bar colors (color any samples with)
            colors = ['gray' if val >= recommended_min else 'lightgrey' for val in num_sites_in_sample.values]
        else:
            colors = 'gray'
            
                
        if plot_type == 'bars':
            #plot a bar graph
            ax.bar(xticks, num_sites_in_sample, width = self.offset*2, color = colors, edgecolor='black')
            ax.set_ylabel('Evidence Size', rotation = 0, ha = 'right', va = 'center')
            
            #add recommended labels
            if include_recommendations:
                if bar_line_colors is None:
                    bar_line_colors = [palette[8], palette[3]]
                elif len(bar_line_colors) != 2:
                    print('Must provide 2 colors in the following order: ideal min, recommended min. Using default colors')
                    bar_line_colors = [palette[8], palette[3]]
                ax.axhline(ideal_min, c = bar_line_colors[0], linestyle = 'dashed', linewidth = 0.8, label =f'Ideal Minimum (n>{ideal_min})')
                ax.axhline(recommended_min, c = bar_line_colors[1], linewidth = 0.8, label = f'Recommended Minimum (n>{recommended_min})')
                ax.legend(bbox_to_anchor = (1, 1))
                    
        elif plot_type == 'dots':
            colors = []
            if include_recommendations:
                if dot_colors is None:
                    dot_colors = [palette[2],palette[8],palette[0]]
                elif len(dot_colors) != 3:
                    print('Must provide 3 colors in the following order: ideal, sufficient, low. Using default colors')
                    palette = sns.color_palette('colorblind')
                    dot_colors = [palette[2],palette[8],palette[0]]
                for size in num_sites_in_sample:
                    if size > ideal_min:
                        colors.append(dot_colors[0])
                    elif size > recommended_min:
                        colors.append(dot_colors[1])
                    else:
                        colors.append(dot_colors[2])
                        
                #create legend
                legend_elements = [Line2D([0],[0], color = 'white',markerfacecolor = dot_colors[0], 
                                          label = f'Ideal Evidence Size (n>{ideal_min})', marker = 'o', markersize = 10),
                                  Line2D([0],[0], color = 'white',markerfacecolor = dot_colors[1], 
                                         label = f'Sufficient Evidence Size (n>{recommended_min})', marker = 'o', markersize = 10),
                                  Line2D([0],[0], color='white',markerfacecolor = dot_colors[2], 
                                         label = f'Low Evidence Size (n<{recommended_min})', marker = 'o', markersize = 10)]
                ax.legend(handles = legend_elements, bbox_to_anchor = (1,1))
            else:
                colors = 'gray'
            ax.scatter(xticks, np.repeat(0.5, len(num_sites_in_sample)),s = num_sites_in_sample*dot_size, c = colors)
            ax.axes.get_yaxis().set_visible(False)

    def setup_figure(self, cluster_samples = False, cluster_kinases = False, include_evidence = False, include_context = False, number_contexts = 1):
        #calculate ratios based on figsize
        desired_cluster_height = 0.3 * cluster_samples
        desired_cluster_width = 0.7 * cluster_kinases
        desired_context_height = 0.3 * number_contexts *include_context 
        desired_evidence_height = 0.5 * include_evidence
        desired_dots_height = self.figsize[1] - desired_cluster_height - desired_context_height - desired_evidence_height
        

        # Define subplot types and their inclusion flags/ratios
        row_components = [
            (cluster_samples, desired_cluster_height),
            (include_context, desired_context_height),
            (True, desired_dots_height),  # main dotplot always included
            (include_evidence, desired_evidence_height)
        ]

        desired_cluster_width = 0.3 * cluster_kinases
        desired_dots_width = self.figsize[0] - desired_cluster_width
        col_components = [
            (cluster_kinases, desired_cluster_width),
            (True, desired_dots_width)  # main dotplot always included
        ]

        # Filter for included components and extract ratios
        height_ratios = [ratio for include, ratio in row_components if include]
        width_ratios = [ratio for include, ratio in col_components if include]


        nrows = len(height_ratios)
        ncols = len(width_ratios)

        fig, axes = plt.subplots(figsize = self.figsize, nrows=nrows, ncols=ncols, gridspec_kw={'height_ratios': height_ratios, 'width_ratios': width_ratios}, sharex = 'col', sharey = 'row')
        fig.subplots_adjust(wspace=0, hspace = 0)

        #if multiple columns and rows, remove unused axes (in top left corner of plot)
        if nrows > 1 and ncols > 1:
            if include_evidence:
                for i in range(nrows - 2):
                    axes[i, 0].axis('off')
                axes[-1,0].axis('off')
            else:
                for i in range(nrows - 1):
                    axes[i, 0].axis('off')
            
            axes[-1, 0].set_xticks([])

        elif ncols > 1:
            axes[0].set_xticks([])


        row = include_context + cluster_samples
        if cluster_kinases and any([cluster_samples, include_context, include_evidence]):
            dots_ax = axes[row, 1]
        elif any([cluster_samples, include_context, include_evidence]):
            dots_ax = axes[row]
        elif cluster_kinases:
            dots_ax = axes[1]
        else:
            dots_ax = axes
        dots_ax.tick_params(axis = 'x', rotation = 90)
        return fig, axes, nrows, ncols, dots_ax


    def make_complete_dotplot(self, kinases_to_plot=None, cluster_samples = False, cluster_kinases = False, sort_kinases_by = None, sort_samples_by = None, binary_evidence = None, context = None, significant_kinases_only = True, show_xtick_labels = True, **kwargs):
        """
        Master function for creating a comprehensive dotplot visualization, which automatically creates any necessary subplots

        Parameters
        ----------
        kinases_to_plot : list or None, optional
            List of kinases to include in the plot. If None, all kinases are included.
        cluster_samples : bool, optional
            Whether to cluster samples in the plot.
        cluster_kinases : bool, optional
            Whether to cluster kinases in the plot.
        significant_kinases_only : bool, optional
            Whether to include only significant kinases in the plot.
        sort_samples_by : str or None, optional
            Kinase Column to sort samples by in the plot based on kinase activities. If cluster_sample=True, this will be ignored.
        sort_kinases_by : str or None, optional
            Sample Column to sort kinases by in the plot based on kinase activities. If cluster_kinases=True, this will be ignored.
        binary_evidence : pd.DataFrame or None
            Binary evidence dataframe from KSTAR analysis. If provided, will calculate the number of sites used as evidence in each sample and plot this.
        context : pd.DataFrame or None, optional
            Context dataframe providing additional sample information for plotting. If provided, must include an 'id_column' for unique sample identifiers and list 'context_columns' for context information.
        show_xtick_labels : bool, optional
            Whether to show x-axis tick labels in the dotplot.
        **kwargs : 
            Additional keyword arguments passed to plotting functions, like matplotlib.pyplot.scatter, DotPlot.context, DotPlot.dotplot, DotPlot.cluster, and DotPlot.evidence_count
        """

        #process kwargs, report any that are not used
        scatterplot_kwargs = helpers.extract_kwonlyargs(plt.scatter, **kwargs)
        #cluster kwargs
        cluster_kwargs = helpers.extract_relevant_kwargs(self.cluster, **kwargs) if cluster_samples or cluster_kinases else {}
        #context kwargs, combine with scatterplot
        context_kwargs = helpers.extract_relevant_kwargs(self.context, **kwargs) if context is not None else {}
        context_kwargs = {**context_kwargs, **scatterplot_kwargs}
        #evidence_count kwargs
        evidence_kwargs = helpers.extract_relevant_kwargs(self.evidence_count, **kwargs) if binary_evidence is not None else {}
        #dotplot kwargs, combine with scatterplot kwargs
        dotplot_kwargs = helpers.extract_relevant_kwargs(self.dotplot, **kwargs)
        dotplot_kwargs = {**dotplot_kwargs, **scatterplot_kwargs}


        #identify any unused kwargs
        used_kwargs = set(cluster_kwargs) | set(context_kwargs) | set(evidence_kwargs) | set(dotplot_kwargs)
        unused_kwargs = {k: v for k, v in kwargs.items() if k not in used_kwargs}
        if unused_kwargs:
            print(f"Warning: The following kwargs were not recognized: {unused_kwargs}")

        

        #check if binary evidence is provided and should be included in plot
        include_evidence = 1 if binary_evidence is not None else 0

        #make sure necessary keyword arguments are set for context
        if context is not None:
            id_column = context_kwargs.get('id_column', None)
            context_columns = context_kwargs.get('context_columns', None)
            if id_column is None or context_columns is None:
                raise ValueError("Both 'id_column' and 'context_columns' must be provided in kwargs when context is specified. id_column will indicate the column with unique sample identifiers, and context_columns will indicate the columns with context information.")
            elif not isinstance(context_columns, list):
                raise ValueError("'context_columns' must be provided as a list of all columns with context information. If only one type of context is provided, it should still be in a list.")
            
            include_context =1
            number_contexts = len(context_columns)
        else:
            include_context = 0
            number_contexts = 0


        if kinases_to_plot:
            kinases_to_drop = [kin for kin in self.values.index if kin not in kinases_to_plot]
            if kinases_to_drop:
                self.drop_kinases(kinases_to_drop)
            #report any kinases that are not recognized/in the activities dataframe
            unrecognized_kinases = [kin for kin in kinases_to_plot if kin not in self.values.index]

            if unrecognized_kinases:
                print(f"Warning: The following kinases are not recognized and will be ignored: {unrecognized_kinases}")

        if significant_kinases_only:
            self.drop_kinases_with_no_significance()
            #if specific kinases were requested but not significant, report
            if kinases_to_plot is not None:
                insignificant_kinases = [kin for kin in kinases_to_plot if kin not in self.values.index and kin not in unrecognized_kinases]
                if insignificant_kinases:
                    print(f"Warning: The following kinases were requested but not significant and will be ignored: {insignificant_kinases}. Set significant_kinases_only=False to include them.")

        if cluster_kinases and self.values.shape[0] <= 1:
            cluster_kinases = False
            print("Warning: Cannot cluster kinases when only one kinase is present.")
        if cluster_samples and self.values.shape[1] <= 1:
            cluster_samples = False
            print("Warning: Cannot cluster samples when only one sample is present.")

        fig, axes, nrows, ncols, dots_ax = self.setup_figure(cluster_samples=cluster_samples, cluster_kinases=cluster_kinases, include_evidence=include_evidence, include_context=include_context, number_contexts=number_contexts)

        #sort by kinase, if provided
        if sort_kinases_by is not None and not cluster_kinases:
            self.values = self.values.sort_values(by=sort_kinases_by)
            self.fpr = self.fpr.loc[self.values.index]
            self.set_index_labels(self.kinase_dict) #update index labels after sorting
            self.set_colors() #update significance colors after sorting
        elif sort_kinases_by is not None and cluster_kinases:
            print(f"Warning: Sorting kinases by a {sort_kinases_by} will not be applied when clustering kinases.")

        #cluster the kinases with hierarchical clustering if desired
        if cluster_kinases:
            row = 0 + include_context + cluster_samples
            col = 0
            if nrows > 1:
                ax = axes[row, col]
            else:
                ax = axes[row]

            self.cluster(orientation = 'left', ax = ax, method='ward', **cluster_kwargs)


        #sort by sample, if provided
        if sort_samples_by is not None and not cluster_samples:
            self.values = self.values.sort_values(by=sort_samples_by, axis=1)
            self.fpr = self.fpr.loc[:, self.values.columns]
            self.set_column_labels(self.x_label_dict) #update column labels after sorting
            self.set_colors() #update significance colors after sorting
        elif sort_samples_by is not None and cluster_samples:
            print(f"Warning: Sorting samples by a {sort_samples_by} will not be applied when clustering samples.")


        if cluster_samples:
            ax = axes[0,1] if cluster_kinases else axes[0]
            self.cluster(orientation = 'top', ax = ax, method='ward', **cluster_kwargs)
            axes[0,0].axis('off')


        if context is not None:
            row = 0 + cluster_samples
            if cluster_kinases:
                ax = axes[row, 1]
            else:
                ax = axes[row]


            self.context(ax = ax, info = context, orientation = 'top', **context_kwargs)

        if binary_evidence is not None:
            row = 1 + include_context + cluster_samples
            ax = axes[row, 1] if cluster_kinases else axes[row]

            self.evidence_count(ax = ax, binary_evidence = binary_evidence, **evidence_kwargs)
            if cluster_kinases:
                axes[nrows-1,0].axis('off')
            ax.tick_params(axis='x', rotation=90)



        self.dotplot(ax=dots_ax, **dotplot_kwargs)
        if not show_xtick_labels:
            dots_ax.set_xticks([])
        




class KSTAR_PDF(fpdf.FPDF):
    """
    Class to generate a PDF report from KSTAR analysis results.
    """
    def __init__(self, activities, fpr, odir, name, binarized_experiment, param_dict):
        super().__init__()
        #detect if values are already log-transformed. If not, transform them
        if np.all((activities >= 0) & (activities <= 1)):
            activities = -np.log10(activities)

        self.activities = activities
        #convert data columns to readable strings by fpdf
        self.activities.columns = [self.filter_unicode_chars(col) for col in self.activities.columns]
        self.fpr = fpr
        self.fpr.columns = [self.filter_unicode_chars(col) for col in self.fpr.columns]
        self.binarized_experiment = binarized_experiment
        self.binarized_experiment.columns = [self.filter_unicode_chars(col) for col in self.binarized_experiment.columns]

        #check if there are special characters in the column names that may cause issues with font, remove them if so

        self.odir = odir
        self.name = name
        self.param_dict = param_dict
        self.phospho_type = param_dict['phospho_type']
        #check if Figures directory exists, if not create it
        if not os.path.exists(os.path.join(self.odir, "FIGURES")):
            os.makedirs(os.path.join(self.odir, "FIGURES"))

    def table(self, data, header = None,column_widths = 40, row_height = 5):
        if header is None:
            header = data.columns
        # Column widths
        self.set_font('Helvetica', 'B', 10)
        
        # Header
        for i, header_item in enumerate(header):
            if isinstance(column_widths, int):
                self.cell(column_widths, row_height, header_item, 1)
            elif isinstance(column_widths, list):
                self.cell(column_widths[i], row_height, header_item, 1)
            else:
                raise ValueError('column_widths must be an integer or a list of integers')

        self.ln(row_height)
        
        # Data
        self.set_font('Helvetica', '', 10)
        for i, row in data.iterrows():
            for i, item in enumerate(row):
                if isinstance(column_widths, int):
                    self.cell(column_widths, row_height, str(item), 1)
                elif isinstance(column_widths, list):
                    self.cell(column_widths[i], row_height, str(item), 1)
                else:
                    raise ValueError('column_widths must be an integer or a list of integers')
            self.ln(row_height)

    def filter_unicode_chars(self, text):
        if not isinstance(text, str):
            text = str(text)
        allowed_chars = set(string.printable)

        #below translator code taken directly from stack overflow post https://stackoverflow.com/questions/59552782/how-to-convert-characters-from-greek-to-english-python)
        greek_alphabet = 'ΑαΒβΓγΔδΕεΖζΗηΘθΙιΚκΛλΜμΝνΞξΟοΠπΡρΣσςΤτΥυΦφΧχΨψΩω'
        latin_alphabet = 'AaBbGgDdEeZzHhJjIiKkLlMmNnXxOoPpRrSssTtUuFfQqYyWw'
        greek2latin = str.maketrans(greek_alphabet, latin_alphabet)
        filtered_text = text.translate(greek2latin)

        #remove any remaining non-allowed characters
        filtered_text = ''.join(c for c in filtered_text if c in allowed_chars)
        return filtered_text

    def mapping_page(self):
        pass


    def summary_page(self):
        self.add_page()
        # Helvetica bold 15
        self.set_font('Helvetica', 'B', 15)
        # Move to the right
        self.cell(80)
        self.cell(20, 10, f'Summary of KSTAR Run ({self.param_dict["phospho_type"]})', 0, 0, 'C')

        #report parameters used
        self.ln(10)
        if self.param_dict is not None:
            self.set_font('Helvetica', 'B', 14)
            self.cell(0, 5, 'Parameters Used:', 0, 1, 'L')
            self.ln(5)
            #iterate through parameters and print them
            for key, value in self.param_dict.items():
                if isinstance(value, (bool, float, int, str)) and key not in ['network_check', 'network_directory', 'pregenerated_experiments_path', 'mann_whitney', 'kinases', 'randomized']:
                    self.set_font('Helvetica', 'B', 10)
                    #set cell width to 
                    self.cell(w=53, h=4, text=f'{self.filter_unicode_chars(key)}: ', border=0, new_x='RIGHT', new_y='TOP', align='R')
                    self.set_font('Helvetica', '', 10)
                    self.cell(w=0, h=4, txt=f'{self.filter_unicode_chars(value)}', border=0, new_x='LMARGIN', new_y='NEXT', align='L')



            self.ln(10)
        
        # Summary statistics
        self.set_font('Helvetica', 'B', 14)
        self.cell(w=0, h=4, txt='Summary Statistics:', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        self.ln(5)
        self.set_font('Helvetica', '', 10)
        num_samples = self.activities.shape[1]
        num_kinases = self.activities.shape[0]
        self.cell(w=0, h=4, txt=f'Number of Samples/Columns Analyzed: {num_samples}', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        self.cell(w=0, h=4, txt=f'Number of Kinases Analyzed: {num_kinases}', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        significant_kinases = (self.fpr < 0.05).any(axis = 1).sum()
        self.cell(w=0, h=4, txt=f'Number of Kinases with significant activity in at least one sample: {significant_kinases}', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        self.ln(10)
        self.top_kinases_table()
        

    def create_dotplot(self):
        print('generating dotplot figure...')
        #calculate number of samples and kinases to size figure appropriately
        num_samples = self.activities.shape[1]
        num_kinases = (self.fpr < 0.05).any(axis = 1).sum()
        if num_kinases*0.2 < 3:
            fig_height = 3
        elif num_kinases*0.2 >= 12:
            fig_height = 12
        else:
            fig_height = 0.2 * num_kinases
        #fig_width = 0.5 * num_samples if 0.5 * num_samples < 8 else 8
        fig_width = 8


        #construct dotplot
        dot = DotPlot(self.activities, self.fpr, legend_title='-log10(p-value)', figsize = (fig_width, fig_height))

        if self.param_dict is not None:
            phospho_type = self.param_dict['phospho_type']
            include_recommendations = True
        else:
            phospho_type = None
            include_recommendations = False

        show_xtick_labels = True if num_samples <=50 else False
        dot.make_complete_dotplot(cluster_kinases = True, binary_evidence = self.binarized_experiment, significant_kinases_only=True, phospho_type=phospho_type, include_recommendations=include_recommendations, show_xtick_labels=show_xtick_labels)

        
        plt.savefig(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_dotplot.png"), bbox_inches='tight', transparent=True, dpi=300)
        plt.close()

    def dotplot_page(self, regenerate_plots = False):
        self.add_page()
        # Helvetica bold 15
        self.set_font('Helvetica', 'B', 15)
        # Move to the right
        self.cell(80)
        self.cell(w=30, h=10, txt=f'KSTAR Dotplot', border=0, new_x='LMARGIN', new_y='NEXT', align='C')
        phospho_type = self.param_dict['phospho_type']

        #add indicators of where to find the dotplot figure and underlying data
        self.set_font('Helvetica', 'I', 8)
                #try including data
        #if embed_files:
        #    phospho_type = self.param_dict['phospho_type']
        #    self.file_attachment_annotation(x=3, y=self.get_y(), w=5, h=4, file_path=os.path.join(self.odir, "FIGURES", f"{self.name}_dotplot.png"))
        self.cell(w=5, h=4, txt=f'The KSTAR dotplot can be found in output directory under "{os.path.join("FIGURES", f"{self.name}_dotplot.png")}".', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        #if embed_files:
        #    phospho_type = self.param_dict['phospho_type']
        #    self.file_attachment_annotation(x=3, y=self.get_y(), w=5, h=4, file_path=os.path.join(self.odir, "RESULTS", f"{self.name}_{phospho_type}_mann_whitney_activities.tsv"))
        self.cell(w=0, h=4, txt=f'Raw activity scores (dot size) in the output directory under "{os.path.join("RESULTS", f"{self.name}_{phospho_type}_mann_whitney_activities.tsv")}".', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        #if embed_files:
        #    phospho_type = self.param_dict['phospho_type']
        #    self.file_attachment_annotation(x=3, y=self.get_y(), w=5, h=4, file_path=os.path.join(self.odir, "RESULTS", f"{self.name}_{phospho_type}_mann_whitney_fpr.tsv"))
        self.cell(w=0, h=4, txt=f'False positive rates (significance) can be found in "{os.path.join("RESULTS", f"{self.name}_{phospho_type}_mann_whitney_fpr.tsv")}".', border=0, new_x='LMARGIN', new_y='NEXT', align='L')

        #create and add dotplot to page
        if not os.path.exists(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_dotplot.png")) or regenerate_plots:
            self.create_dotplot()
        self.image(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_dotplot.png"), x = 0, y = self.get_y()+10, w = 200)

        #add link to kstar plotting tool
        self.ln(10)
        self.set_font('Helvetica', 'B', 12)
        self.set_y(-30)
        self.cell(w=0, h=5, txt='Want to customize the dotplot? Visit the KSTAR plotting tool linked here:', border=0, new_x='LMARGIN', new_y='NEXT', align='L', link='https://proteomescout.research.virginia.edu/kstar/')

    def evidence_count_plot(self, data_columns):
        #get evidence counts
        data_columns = self.activities.columns.tolist()
        num_samples = len(data_columns)
        num_sizes = self.binarized_experiment[data_columns].sum()
        #remove 'data:' from column names for plotting
        num_sizes.index = [col.replace('data:','') for col in num_sizes.index]

        #make figure
        fig, ax = plt.subplots(figsize = (7,3.5))
        ax.barh(y = num_sizes.index, width = num_sizes.values, color = 'gray', edgecolor='black', height = 1, lw = 0.5)
        ax.set_xlabel('Number of Sites Used as Evidence')
        ax.set_ylabel('Column')

        #adjust y-axis label size based on number of samples
        if num_samples > 20 and num_samples <= 40:
            ax.tick_params(axis = 'y', labelsize=6)
        elif num_samples > 40:
            ax.set_yticks([])
        
        if num_samples <= 35:
            #annotate bars with counts
            for i, v in enumerate(num_sizes.values):
                ax.text(v + 0.5, i, str(int(v)), color='black', va='center', ha = 'left', fontsize=7)

        #add vertical line for median evidence size
        median_size = np.median(num_sizes.values)
        ax.axvline(median_size, color = 'red', linestyle = 'dashed', linewidth = 1, label = f'Median Evidence Size ({int(median_size)})')
        ax.legend(loc = (0, 1.05))

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        #save and close
        plt.savefig(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_count.png"), bbox_inches='tight', transparent=True, dpi=300)
        plt.close()

    def evidence_overlap_plot(self, data_columns):
        fig, ax = plt.subplots(figsize = (7,3.5))
        if len(data_columns) <= 20:
            annot = True
            xticklabels = 1
            yticklabels = 1
            xticksize = 9
            yticksize = 9
        elif len(data_columns) <= 40:
            annot = False
            xticklabels = 1
            yticklabels = 1
            xticksize = 6
            yticksize = 6
        else:
            annot = False
            xticklabels = False
            yticklabels = False
            xticksize = 6
            yticksize = 6

        if self.activities.shape[1] > 2:
            ax = plot_jaci_between_samples(self.binarized_experiment, data_columns, ax = ax, annot_kws={"size": 6}, cluster = True, linewidths = 0.1, linecolor = 'black', annot=annot, xticklabels = xticklabels, yticklabels = yticklabels)
        else:
            ax = plot_jaci_between_samples(self.binarized_experiment, data_columns, ax = ax, annot_kws={"size": 6}, cluster = False, linewidths = 0.1, linecolor = 'black', annot=annot, xticklabels = xticklabels, yticklabels = yticklabels)

        ax.tick_params(axis = 'x', labelsize=xticksize)
        ax.tick_params(axis = 'y', labelsize=yticksize)
        plt.savefig(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_overlap.png"), bbox_inches='tight', transparent=True, dpi=300)
        plt.close()

    def evidence_page(self, regenerate_plots = False):
        data_columns = self.activities.columns.tolist()
        
        if not os.path.exists(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_count.png")) or regenerate_plots:
            self.evidence_count_plot(data_columns)
        
        if not os.path.exists(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_overlap.png")) or regenerate_plots:
            self.evidence_overlap_plot(data_columns)

        self.add_page()
        self.set_font('Helvetica', 'B', 15)
        self.cell(80)
        self.cell(w=30, h=10, txt='Evidence Characteristics', border=0, new_x='LMARGIN', new_y='NEXT', align='C')

        #add indicators of where to find the dotplot figure and underlying data
        self.set_font('Helvetica', 'I', 8)
        self.cell(w=0, h=4, txt=f'The evidence barplot can be found at {os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_count.png")}.', border=0, new_x='LMARGIN', new_y='NEXT', align='L')
        self.cell(w=0, h=4, txt=f'The sites used as evidence can be found at {os.path.join(self.odir, "RESULTS", f"{self.name}_{self.phospho_type}_binarized_experiment.tsv")}.', border=0, new_x='LMARGIN', new_y='NEXT', align='L')

        #plot evidence sizes as a barplot
        self.ln(10)
        self.set_font('Helvetica', 'B', 12)
        self.cell(0, 10, 'Evidence Sizes Across Samples:', 0, 1, 'L')
        self.image(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_count.png"), x = 10, y = self.get_y(), w = 180)

        #plot evidence overlap as a heatmap
        self.ln(115)
        self.set_font('Helvetica', 'B', 12)
        self.cell(0, 10, 'Evidence Overlap Between Samples (By Jaccard Index):', 0, 1, 'L')
        self.image(os.path.join(self.odir, "FIGURES", f"{self.name}_{self.phospho_type}_evidence_overlap.png"), x = 10, y = self.get_y(), w = 180)



    def top_kinases_table(self):
        #construct table of top kinases per sample
        data_columns = self.activities.columns
        num_sig_list = []
        top5_list = []
        #iterate through columns to get statistics
        for col in data_columns:
            significant_kinases = self.fpr[self.fpr[col] < 0.05].index.tolist()
            #get number of significant kinases
            num_significant = len(significant_kinases)
            num_sig_list.append(num_significant)
            #get top 5 kinases by activity
            top5 = self.activities.loc[significant_kinases, col].sort_values(ascending = False).head(5).index.tolist()
            top5_list.append(', '.join(top5))

        #trim 'data:' from column names
        data_columns = [col.replace('data:','') for col in data_columns]
        #convert common greek letters to 
        top_kinases = pd.DataFrame({'Sample':data_columns, '# of Sig. Kinases':num_sig_list, 'Top 5 Kinases':top5_list})

        self.table(top_kinases, column_widths = [50, 30, 100], row_height=6)

        # Page footer
    def footer(self):
        """
        Override the footer method to add a page number at the bottom center of each page.
        """
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Helvetica italic 8
        self.set_font('Helvetica', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')

    def generate(self, regenerate_plots = False):
        #make figure directory if it does not exist
        if not os.path.exists(os.path.join(self.odir, 'FIGURES')):
            os.makedirs(os.path.join(self.odir, 'FIGURES'))

        self.summary_page()
        self.dotplot_page(regenerate_plots=regenerate_plots)
        #self.top_kinases_page()
        self.evidence_page(regenerate_plots=regenerate_plots)



        fname = os.path.join(self.odir, f'{self.name}_{self.phospho_type}_summary.pdf')
        self.output(fname)

        