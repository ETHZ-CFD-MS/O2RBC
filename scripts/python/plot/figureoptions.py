"""
Module containing the class FigureOptions.
"""

import os
import subprocess
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


class FigureOptions:
    """
    Wrapper for figure options related to size, fonts and image quality.
    Encapsulates saving of figures, logfile output and whether to show and save plots.
    Can modify figure (e.g., transform to black/white).
    """
    
    # Useful link:
    # http://wiki.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    inches_per_pt = 1.0/72.27
    inches_per_cm = 1.0/2.54
    golden_mean   = (np.sqrt(5)-1.0)/2.0  # nice looking ratio

    # TODO: write a version with a classmethod to construct from a parser
    # def __init__(self, **kwargs):
    #     self.unit          = 'inch'
    #     self.fig_width     = 3.5
    #     self.height_ratio  = self.golden_mean  # height/width ratio
    #     self.dpi           = 150.0        # dot per inches for saving figure
    #
    #     self.fileFormats   = ['png']
    #     self.textSize      = 'medium'
    #     self.hasLog        = False
    #     self.show          = True
    #     self.saveFile      = True
    #     self.noLegend      = False
    #     self.blackWhite    = False
    #
    #     self.labelpad      = 3
    #
    #     for key in kwargs:
    #         if hasattr(self, key):
    #             setattr(self, key, kwargs.get(key))
    #         else:
    #             warnings.warn('Argument "{:s}" not recognized.'.format(key))

    def __init__(self, parser):
        self.unit          = 'inch'
        self.fig_width     = 3.5
        self.height_ratio  = self.golden_mean  # height/width ratio
        self.dpi           = 150.0        # dot per inches for saving figure

        self.fileFormats   = ['png']
        self.textSize      = 'medium'
        self.hasLog        = False
        self.show          = True
        self.saveFile      = True
        self.noLegend      = False
        self.blackWhite    = False

        self.labelpad      = 3

        self.addOptions(parser)

    def addOptions(self, parser):
        group = parser.add_argument_group('figureOptions')
        group.add_argument('--width', help='Figure width (default unit is inches)', 
                           type=float, default=self.fig_width)
        group.add_argument('--figUnit', help='Figure length unit (inch, cm, pt)', 
                           default=self.unit)
        group.add_argument('--ratio', help='Ratio height/width', 
                           type=float, default=self.golden_mean)
        group.add_argument('--dpi', help='Dot per inches', 
                           type=float, default=self.dpi)
        group.add_argument('--log', help='Write a log file for the figure', 
                           action='store_true')
        group.add_argument('--show', help='Show the plot after saving',
                           action='store_true')
        group.add_argument('--noSave', help='Do not save any file (overrides --log)',
                           action='store_true')
        group.add_argument('--formats', help='Formats in which figure is saved', 
                           nargs='+', default=self.fileFormats)
        group.add_argument('--suffix', help='Suffix for the plot name', default='')
        group.add_argument('--bw',  help='Plot in black and white', 
                           action='store_true')
        group.add_argument('--noLegend',  help='Removes the plot legend',
                           action='store_true')

    def parseOptions(self,args):
        self.unit = args.figUnit        
        self.fig_width = args.width
        # convert length unit to inches
        if self.unit == 'pt':
            self.fig_width = self.fig_width*self.inches_per_pt
        elif self.unit == 'cm':
            self.fig_width = self.fig_width*self.inches_per_cm
        elif self.unit != 'inch':
            raise ValueError('Invalid unit %s' % self.unit)

        self.height_ratio = args.ratio
        self.dpi          = args.dpi
        self.fileFormats  = args.formats
        self.hasLog       = args.log
        self.saveFile     = not args.noSave 
        self.show         = args.show
        self.suffix       = args.suffix
        self.blackWhite   = args.bw
        self.noLegend     = args.noLegend

        self.applyOptions()

    def applyOptions(self):
        """Updates rcParams with required options.

        This should be called before the desired figure is created.
        """
        fig_height = self.fig_width * self.height_ratio
        fig_size = (self.fig_width, fig_height)
        print 'Figure size = ', fig_size
        if self.textSize == 'medium':
            labelsize = 10
            textsize = 10
            legendsize = 10
            ticklabelsize = 10
        elif self.textSize == 'small':
            labelsize = 8
            textsize = 8
            legendsize = 6
            ticklabelsize = 6

        params = {'axes.labelsize':  labelsize,
                  'font.size':       textsize,
                  'legend.fontsize': legendsize,
                  'xtick.labelsize': ticklabelsize,
                  'ytick.labelsize': ticklabelsize,
                  'axes.labelpad':   0,
                  'text.usetex':     True,
                  'text.latex.preamble': [r'\usepackage{txfonts}', # for upright greek letters
                                          r'\usepackage{lmodern}'],
                  'font.family':     'serif',
                  # 'font.family' :  'lmodern',
                  'font.weight':     'normal',
                  'font.style':      'normal',
                  'figure.figsize':  fig_size,
                  'figure.dpi':      self.dpi}
        plt.rcParams.update(params)

    def adjustAxes(self,pad=0.4):
        for ax in plt.gcf().axes:
            ax.xaxis.labelpad = self.labelpad
            ax.yaxis.labelpad = self.labelpad
            plt.tight_layout(pad=pad)

    def setGrid(self):
        xgridlines = plt.getp(plt.gca(), 'xgridlines')
        ygridlines = plt.getp(plt.gca(), 'ygridlines')
        plt.setp(xgridlines, 'dashes',(1,3), 'linewidth', 0.25)
        plt.setp(ygridlines, 'dashes',(1,3), 'linewidth', 0.25)

    def saveFig(self, plotName):
        """Save the current figure to a given file with given formats
        
        Args:
            plotName: Name of the plot, without extension
        """
        self.adjustAxes()
        if self.suffix:
            plotName += self.suffix
        if self.blackWhite:
            fig = plt.gcf()
            plotName += '_bw'

            # Taken from http://www.david-zwicker.de/code/movie_making.py
            def get_filter(name):
                return lambda x: hasattr(x, 'set_%s'%name) and hasattr(x, 'get_%s'%name)   
            # Taken from http://stackoverflow.com/questions/12201577/convert-rgb-image-to-grayscale-in-python
            def rgb2gray(rgb):
                gray = rgb[3]*np.dot(rgb[:3], [0.299, 0.587, 0.114])
                return '%g' % gray

            def is_gray(color):
                """
                Return True if the given color string is black or a gray level.
                """
                if color == 'k':
                    return True
                else:
                    try:
                        float(color)
                        return True
                    except (ValueError, TypeError):
                        pass
                    return False

            # Markers have an edgecolor property
            for o in fig.findobj(get_filter('edgecolor')):
                col = o.get_edgecolor()
                if len(col) == 4:
                    o.set_edgecolor(rgb2gray(col))
                else:
                    try:
                        o.set_edgecolor(rgb2gray(col[0]))
                    except IndexError:
                        print "Could not convert color ", col, " to gray"

            # For facecolor, set color to the grayscale value that corresponds to the facecolor.
            for o in fig.findobj(get_filter('facecolor')):
                col = o.get_facecolor()
                if len(col) == 4:
                    o.set_facecolor(rgb2gray(col))

            # for o in fig.findobj(get_filter('markeredgecolor')):
            #     col = o.get_markeredgecolor()
            #     if len(col) == 4:
            #         o.set_markeredgecolor(rgb2gray(col))

            for o in fig.findobj(matplotlib.lines.Line2D):
                if not is_gray(o.get_color()):
                    o.set_color('k')
               
            for o in fig.findobj(matplotlib.text.Text):
                if not is_gray(o.get_color()):
                    o.set_color('k')

        if self.saveFile:
            for fmt in self.fileFormats:
                fileName = '%s.%s' % (plotName, fmt)
                if fmt != 'tif' and fmt != 'tiff':
                    plt.savefig(fileName, format=fmt, dpi=self.dpi)
                    print 'Saved plot %s with %i dpi' % (fileName, self.dpi)
                else:
                    from PIL import Image
                    import cStringIO
                    # png1 = cStringIO.StringIO()
                    # print png1
                    tmp_file_name = '.tmp.png'
                    plt.savefig(tmp_file_name, format='png', dpi=self.dpi)
                    png2 = Image.open(tmp_file_name)
                    png2.save(fileName)
                    print 'Saved plot %s with %i dpi' % (fileName, self.dpi)
                    # png1.close()


            if self.hasLog:
                self.writeLog(plotName)

        if self.show:
            plt.show()

    def writeLog(self, plotName):
        """Writes a log file with the figure options and date information.
        
        Args:
            plotName: Name of the plot, without extension
        """
        logName = '%s.txt' % plotName

        with open(logName, 'w') as f:
            f.write('This file is automatically generated by figureoptions.py.\n')
            f.write('\n')
            # write time
            f.write('Figure generated on %s\n' % time.strftime('%c'))
            f.write('Figure name: %s\n' % logName)
            f.write('\n')

            # write figure options
            f.write("Figure width : %g %s\n" % (self.fig_width, self.unit))
            f.write("Figure height: %g %s\n" % 
                                (self.fig_width*self.height_ratio, self.unit))
            f.write("Height ratio : %g\n" % self.height_ratio)
            f.write("dpi          : %g dpi\n" % self.dpi)
            f.write("Text size    : %s\n" % self.textSize)
            f.write("Black-white  : %s\n" % self.blackWhite)
            f.write("No legend    : %s\n" % self.noLegend)
            f.write('\n')
            f.write("Label pad    : %s\n" % self.labelpad)
            f.write('\n')

            # write current path
            f.write('Path: %s\n' % os.getcwd())
            f.write('\n')

            # write revision of git code
            f.write('Last code commit:\n')
            proc = subprocess.Popen('git log -1 | head -n 3', 
                                    stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            f.write(out)

        print 'Wrote log file %s' % logName


