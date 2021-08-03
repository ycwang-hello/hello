一些杂七杂八的程序。

# Solid State Physics

## neighbors

程序是neighbors.py，计算面心立方和体心立方的最近邻、第二近邻、第三近邻、……的格点个数及其坐标，结果在bcc_data.csv和fcc_data.csv中。画了一些示意图，但是立体的有点难想像，看这些图可能看不懂。



# Quantum Mechanics

## Zeeman_plot

Prerequisites: Python 3; matplotlib; seaborn; numpy.

Usage: input L1, S1, J1, L2, S2, J2 of the two energy levels.
Examples: 

```Python
#plot(s1, l1, j1, s2, l2, j2, name1='state 1 name', name2='state 2 name')
plt.figure(); plot(0,2,2,0,1,1, name1=r'$\mathrm{6s6d{}^1D_2}$', name2=r'$\mathrm{6s6p{}^1P_1}$')
plt.figure(); plot(1,0,1,1,1,2, name1=r'$\mathrm{6s7s{}^3S_1}$', name2=r'$\mathrm{6s6p{}^3P_2}$')
```

Outputs: green.png; yellow.png.

# Astrometry

## SR Doppler Aberration

在狭义相对论框架下演示多普勒效应和光行差。