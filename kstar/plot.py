import numpy as np
import pandas as pd
from enum import Enum


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
                 labelmap = None,
                 legend_title = 'p-value', size_number = 5, size_color = 'gray', 
                 color_title = 'Significant', markersize = 10, 
                 legend_distance = 1.0, figsize = (20,4), title = None,
                 xlabel = True, ylabel = True, x_label_dict = None, kinase_dict = None):


        self.values = values.copy()
        self.fpr = fpr.copy()
        #make sure that fpr dataframe has the same index as values dataframe. If not, reindex
        self.fpr = self.fpr.loc[self.values.index,self.values.columns]
        self.alpha = alpha
        #create binary dataframe that indicates significance based on provided fpr cutoff.
        if inclusive_alpha:
            self.significance = (self.fpr <= alpha) * 1
        else:
            self.significance = (self.fpr < alpha) * 1
        #Assign either fpr or significance to colors dataframe based on 
        self.binary_sig = binary_sig
        if binary_sig:
            self.colors = self.significance
            if labelmap is None:
                if inclusive_alpha:
                    self.labelmap = {0: 'FPR > %0.2f'%(alpha), 1:'FPR <= %0.2f'%(alpha)}
                else: 
                    self.labelmap = {0: 'FPR >= %0.2f'%(alpha), 1:'FPR < %0.2f'%(alpha)}

        else:
            self.colors = self.fpr
      
        
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

        self.columns = self.set_column_labels(values, x_label_dict)
        self.index = self.set_index_labels(values, kinase_dict)
    
    def set_values(self, values):
        self.values = values

    def set_colors(self, colors):
        self.colors = colors

    def set_column_labels(self, values, x_label_dict):
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
            
    def set_index_labels(self, values, kinase_dict):
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
        

    def dotplot(self, ax = None, orientation = 'left', size_legend = True, color_legend = True, max_size = None):
        """
        Generates the dotplot plot, where size is determined by values dataframe and color is determined by significant dataframe
        
        Parameters
        -----------
        ax : matplotlib Axes instance, optional
            axes dotplot will be plotted on. If None then new plot generated
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
        scatter = ax.scatter(x, y, c=c, s=s)
        
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
        self.set_column_labels(self.values, self.x_label_dict)
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
            self.set_index_labels(self.values, self.kinase_dict)

        
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
            self.set_column_labels(self.values, self.x_label_dict)

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
        
        
    def context(self, ax, info, id_column, context_columns, dotsize = 200, markersize = 20, orientation = 'left', color_palette='colorblind', margin = 0.2, make_legend = True):
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
            index = list(self.values.index)
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
            ax.scatter(x = melted['variable'], y = melted['var'],c = melted['value'].map(color_map), s = dotsize)
            ax.tick_params(axis = 'x', rotation = 90)
            ax.axes.get_yaxis().set_visible(False)
            ax.margins(margin, 0.05)
        elif orientation in ['top', 'bottom']:
            ax.scatter(x = melted['var'], y = melted['variable'], c = melted['value'].map(color_map), s = dotsize)
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
                
    def evidence_count(self, ax, binary_evidence, plot_type = 'bars', phospho_type = None, dot_size = 1, include_recommendations = True,
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
                print('Please indicate which phospho_type this plot is for to get appropriate recommendations')
                
            #calculate bar colors (color any samples with)
            colors = ['gray' if val >= recommended_min else 'lightgrey' for val in num_sites_in_sample.values]
        else:
            colors = 'gray'
            
                
        if plot_type == 'bars':
            #plot a bar graph
            ax.bar(xticks, num_sites_in_sample, width = self.offset*2, color = colors)
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
        
        
        
      
        