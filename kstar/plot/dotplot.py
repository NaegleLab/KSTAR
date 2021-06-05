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
    """
    
    
    def __init__(self, values, colors, dotsize = 20, 
                 colormap={0: '#6b838f', 1: '#FF3300'}, 
                 labelmap = {0 : 'Not Significant', 1 : 'Significant'},
                 facecolor = 'white',
                 legend_title = 'p-value', size_number = 5, size_color = 'gray', 
                 color_title = 'Significant', markersize = 10, 
                 legend_distance = 1.0, figsize = (20,4), title = None,
                 xlabel = True, ylabel = True, x_label_dict = None, kinase_dict = None):
        """
        Parameters
        ----------
        values: pandas DataFrame instance
            values to plot 
        colors : pandas DataFrame instance 
            color each value should be plotted as
        dotsize : float, optional
            multiplier to use for scaling size of dots
        colormap : dict, optional
            maps color values to actual color to use in plotting
            default : {0 : '#CFD8DC', 1 : '#FF3300'}, 
        labelmap = 
            maps labels of colors
            default : {0 : 'Not Significant', 1 : 'Significant'},
        facecolor : color, optional
            Background color of dotplot
            default : '#455A64',
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
        """

        self.values = values
        self.colors = colors
        self.figsize =  figsize
        
        self.title = title
        self.xlabel = xlabel
        self.ylabel= ylabel

        self.colormap = colormap
        self.labelmap = labelmap
        
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
                self.x_label_dict[col] = col.strip('data:')
            self.column_labels = [x.strip('data:') for x in self.column_labels]

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
        melt['var'] = melt.apply(lambda row : columns.index(row[1]) * self.multiplier + self.offset, axis = 1)
        
        melt_color = self.colors.melt(id_vars = 'row_index')
        melt_color['var'] = melt_color.apply(lambda row : columns.index(row[1]) * self.multiplier + self.offset, axis = 1)
        self.colors.drop(columns = ['row_index'], inplace = True)

        # Plot Data
        x = melt['var']
        y = melt['row_index']
        
        s = melt.value * self.dotsize
        
        #check to see if more than 2 values are given (fprs). Otherwise get color based on binary significance
        if len(np.unique(self.colors)) > 2:
            cmap = LinearSegmentedColormap.from_list("sig_cmap", list(zip([0,1], [self.colormap[0], self.colormap[1]])))
            norm = Normalize(vmin=0, vmax=2, clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
            #replace 0 with 0.01 to avoid log10 errors, transform the fprs with a log transform
            melt_color.replace(0, 0.01, inplace=True)
            melt_color.value = -np.log10(melt_color.value)
            #get color for each datapoint based on fpr value
            melt_color['color'] = [mapper.to_rgba(l) for l in melt_color.value]
        else:  
            #get color for each datapoint based on significance
            melt_color['color'] = [self.colormap.get(l,'black') for l in melt_color.value]
            
        c = melt_color['color']
        scatter = ax.scatter(x, y, c=c, s=s)
        
        # Add Color Legend
        if color_legend:
            if len(np.unique(self.colors)) > 2:
                #choose which values to show in the legend
                legend_vals = [1, 0.5, 0.05, 0.01]
                legend_color = [mapper.to_rgba(-np.log10(val)) for val in legend_vals]
                #create the legend 
                color_legend = []
                for i in range(len(legend_vals)):
                    color_legend.append(Line2D([0], [0], marker='o', color='w', label=str(legend_vals[i]),
                                markerfacecolor= legend_color[i], markersize=self.markersize))
                legend1 = ax.legend(handles=color_legend, loc=f'upper {orientation}', bbox_to_anchor=(self.legend_distance,1), title = self.color_title)  
                ax.add_artist(legend1) 
            else:
                #create the legend
                color_legend = []
                for color_key in self.colormap.keys():
                    color_legend.append(
                        Line2D([0], [0], marker='o', color='w', label=self.labelmap[color_key],
                                markerfacecolor= self.colormap[color_key], markersize=self.markersize),
                    )     
                legend1 = ax.legend(handles=color_legend, loc=f'upper {orientation}', bbox_to_anchor=(self.legend_distance,1), title = self.color_title)  
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
        
        # Add Additional Plotting Information
        ax.tick_params(axis = 'x', rotation = 90)
        ax.yaxis.set_ticks(np.arange(len(self.values)) * self.multiplier + self.offset)
        ax.xaxis.set_ticks(np.arange(len(columns)) * self.multiplier + self.offset)
        
        # set column labels in case values has changed
        self.set_column_labels(self.values, self.x_label_dict)
        ax.set_xticklabels(self.column_labels)
        ax.set_yticklabels(self.index_labels)
        #adjust yscale so that data is always equally spaced
        ax.set_ylim([0,len(self.values)*self.multiplier])
        ax.set_xlim([0,len(columns)*self.multiplier])
        
        if not self.xlabel:
            ax.axes.xaxis.set_visible(False)
        if not self.ylabel:
            ax.axes.yaxis.set_visible(False)
        return ax 
    
    def cluster(self, ax, method='single', metric='euclidean', orientation = 'top'):
        """
        Performs hierarchical clustering on data and plots result to provided Axes. 
        result and significant dataframes are ordered according to clustering
        
        Parameters
        ---------
        ax : matplotlib Axes instance
            Axes to plot dendogram to
        
        method : str, optional
            The linkage algorithm to use.
        metric : str or function, optional
            The distance metric to use in the case that y is a collection of observation vectors; 
            ignored otherwise. See the pdist function for a list of valid distance metrics. A custom distance function can also be used.
        
        orientation : str, optional
            The direction to plot the dendrogram, which can be any of the following strings:
            'top'
                Plots the root at the top, and plot descendent links going downwards. (default).
            'bottom'
                Plots the root at the bottom, and plot descendent links going upwards.
            'left'
                Plots the root at the left, and plot descendent links going right.
            'right'
                Plots the root at the right, and plot descendent links going left.
        """
        if orientation in ['left', 'right']:
            row_linkage = linkage(self.values, method = method, metric = metric)
            den_row = dendrogram(row_linkage, 
                        ax = ax, 
                        orientation = orientation, 
                        labels = list(self.values.index), 
                        color_threshold = -np.inf, 
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
                            color_threshold = -np.inf, 
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

        kinase_list = self.colors[self.colors.sum(axis=1) ==0].index.values
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

        self.values.drop(index=kinase_list, inplace=True)
        self.colors.drop(index = kinase_list, inplace=True)
        
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
            
        melted = info[[id_column] + context_columns].melt(id_vars=id_column)
        #weird issue with melt function here, where for one datset it provides the context column names in 0 column rather than 'variable'. Rename for now.
        if 0 in melted.columns:
            melted.rename(columns = {0: 'variable'}, inplace = True)
        melted['var'] = melted.apply(lambda row : index.index(row[0]) * self.multiplier + self.offset, axis = 1)
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
                        
                ax.add_artist(leg)
    
                running_total += len(ids) + 1

      
        