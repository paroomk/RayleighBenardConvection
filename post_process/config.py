#! /usr/bin/env python
# -*- coding: utf-8 -*-

# SET UP PLOT PARAMETERS FOR ALL PLOTS
pars = {  'backend'             : 'pdf',
          'text.usetex'         : True,
          'text.fontsize'       : 16,
          'xtick.labelsize'     : 16,
          'xtick.major.pad'     : 10,
          'xtick.major.size'    : 10,
          'xtick.minor.size'    : 5,
          'xtick.major.width'   : 2.5,
          'xtick.minor.width'   : 1.25,
          'ytick.labelsize'     : 16,
          'ytick.major.pad'     : 8,
          'ytick.major.size'    : 10,
          'ytick.minor.size'    : 5,
          'ytick.major.width'   : 2.5,
          'ytick.minor.width'   : 1.25,
          'legend.borderpad'    : 0.25,     # empty space around the legend box
          'legend.labelspacing' : 0.5,
          'legend.borderaxespad': 0.75,
          'legend.fontsize'     : 16,
          'legend.fancybox'     : True,
          'legend.shadow'       : False,
          'legend.frameon'      : True,
          'lines.markersize'    : 7,
          'lines.linewidth'     : 1,
          'font.family'         : 'serif',
          'font.sans-serif'     : 'Helvetica',
          'font.serif'          : 'cmr10',
          'font.size'           : 16,
          'axes.facecolor'      : 'white',   # axes background color
          'axes.edgecolor'      : 'black',   # axes edge color
          'axes.linewidth'      : 1.0,     # edge linewidth
          'axes.labelweight'    : 'normal',  # weight of the x and y labels
          'axes.labelcolor'     : 'black',
          'axes.labelsize'      : 16,
          'axes.axisbelow'      : 'True',
          'axes.formatter.limits': (-4, 4), 
          #'axes.color_cycle'    : ['#95d0fc', '#bf77f6', '#ad8150', '#840000', '#feb308', '#f7022a',
          #                         '#3f9b0b', '#c0fb2d', '#069af3', '#ca6641', '#caa0ff', '#016795',
          #                         '#017b92', '#d3494e', '#f075e6', '#3c9992', '#ff964f', '#6258c4',
          #                         '#6ecb3c', '#aeff6e', '#a442a0', '#0e87cc', '#71aa34', '#fedf08']
          'axes.color_cycle'    : ['#bf5700', '#333f48', '#005f86', '#43695b', '#f2a900', '#382f2d',
                                   '#3f9b0b', '#c0fb2d', '#069af3', '#ca6641', '#caa0ff', '#016795',
                                   '#017b92', '#d3494e', '#f075e6', '#3c9992', '#ff964f', '#6258c4',
                                   '#6ecb3c', '#aeff6e', '#a442a0', '#0e87cc', '#71aa34', '#fedf08']
          }
################################################################################
# COLORS
red            = 'red'
green          = 'green'
blue           = 'blue'
cyan           = 'cyan'
yellow         = 'yellow'
magenta        = 'magenta'
black          = 'black'
orange         = '#CC3300'
purple         = '#700880'
sienna         = '#A0522D'
ForestGreen    = '#228B22'
GoldenRod      = '#DAA520'
CornflowerBlue = '#6495ED'
Crimson        = '#DC143C'
DeepPink       = '#FF1493'
DodgerBlue     = '#1E90FF'
SteelBlue      = '#4682B4'
SkyBlue        = '#87CEEB'
AquaMarine     = '#7FFFD4'
NeonGreen      = '#7FFF00'
DarkMagenta    = '#8B008B'
DarkOrange     = '#FF8C00'
Chocolate      = '#D2691E'
SlateGray      = '#708090'
DarkSlateGray  = '#2F4F4F'
SaddleBrown    = '#8B4513'

# LINESTYLES
solid   = 'solid'
dashed  = 'dashed'
dashdot = 'dash_dot'
dotted  = 'dotted'

# MARKER SHAPES
none = ' '
circle = 'o'
diamond = 'D'
square = 's'
hexagon = 'H'
triright = '>'
trileft = '<'
ex = 'x'
triup = '^'
pentagon = 'p'
star = '*'
clubs = r'$\clubsuit$'
bowtie = r'$\bowtie$'












