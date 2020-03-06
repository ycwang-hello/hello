# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 12:53:57 2019

@author: Yuchen Wang
"""

import matplotlib.pyplot as plt
import numpy as np

import seaborn
seaborn.set_context('paper', font_scale=1.3)

def cal_g(L,S,J):
    g = 1+(J*(J+1)-L*(L+1)+S*(S+1))/(2*J*(J+1))
#    gf = gj*(F*(F+1)+J*(J+1)-I*(I+1))/(2*F*(F+1))
    return g

def plot(s1, l1, j1, s2, l2, j2, name1='', name2='', nl1 = '', nl2 = '', detail=True, spectrum=True):
    g1 = cal_g(l1, s1, j1)
    g2 = cal_g(l2, s2, j2)
    m1s = np.arange(-j1, j1+1)
    m2s = np.arange(-j2, j2+1)
    de1s = g1*m1s
    de2s = g2*m2s
    DE = (de2s[-1]-de1s[0])
    
#    plt.subplot(211)
    #plot 跃迁 
    n = 1
    dx = 0.7
    d_nu_s=[]
    polars = []
    for dm in [1,0,-1]:
        polar = '$\\sigma^-$' if dm == 1 else '$\\sigma^+$' if dm == -1 else '$\\pi$'
        for de1, m1 in zip(de1s, m1s):
            dest = np.argwhere(m2s==m1+dm)
            if not dest.size == 0:
                de2 = de2s[dest][0,0]
                length = 2*DE+de1-de2
#                plt.arrow(n*dx, DE+de1, 0, -length, color='b', width=length*0.0005, head_width = length*0.007, head_length=length*0.07, length_includes_head=True)
                plt.arrow(n*dx, DE+de1, 0, -length, color='b', width=0.005, head_width = 0.08, head_length=DE*0.15, length_includes_head=True)
                plt.text(n*dx, -DE+de2s[0]-g2, polar, horizontalalignment='center')
                plt.text(n*dx, -DE+de2s[0]-2*g2, np.round(de1-de2, 3), horizontalalignment='center')
                d_nu_s.append(de1-de2)
                polars.append(polar)
                
                n+=1
    
    #polar and d nu
    plt.text(-0.3, -DE+de2s[0]-g2, 'pol.', horizontalalignment='center')
    plt.text(-0.3, -DE+de2s[0]-2*g2, '$\\Delta\\tilde{\\nu}/\\tilde{L}$', horizontalalignment='center')
    
    #m and mg
    moreDx = 1.5
    left=-5
    xm = (n)*dx+moreDx+0.5
    xmg = xm+0.8
    plt.text(xm, DE+de1s[-1]+g1, '$m$', verticalalignment='center', horizontalalignment='center')
    plt.text(xmg, DE+de1s[-1]+g1, '$mg$', verticalalignment='center', horizontalalignment='center')
    plt.xlim(left, xmg+0.5)   
    low = -DE+de2s[0]-2*g2-0.5
    high = DE+de1s[-1]+0.5
    plt.ylim(low, high)         
                                
    #plot energy levels
    plt.hlines(np.concatenate((DE+de1s, -DE+de2s)), 0, n*dx+moreDx)
    plt.hlines([-DE, DE], left, -1)
    for de1, m1 in zip(de1s, m1s):
        plt.plot([-1, 0], [DE, DE+de1], '--',color='b')
        plt.text(xm, DE+de1, m1, verticalalignment='center', horizontalalignment='center')
        plt.text(xmg, DE+de1, np.round(de1, 3), verticalalignment='center', horizontalalignment='center')
    for de2, m2 in zip(de2s, m2s):
        plt.plot([-1, 0], [-DE, -DE+de2], '--',color='b')
        plt.text(xm, -DE+de2, m2, verticalalignment='center', horizontalalignment='center')
        plt.text(xmg, -DE+de2, np.round(de2, 3), verticalalignment='center', horizontalalignment='center')
    
    #original energy level names
    plt.text((left-1)/2, DE+0.3, name1, horizontalalignment='center')
    plt.text((left-1)/2, -DE+0.3, name2, horizontalalignment='center')    
    plt.text((left-1)/2, DE-0.3, '$S={}, L={}, J={}, g={}$'.format(s1, l1, j1, np.round(g1, 3)), horizontalalignment='center',verticalalignment='top')    
    plt.text((left-1)/2, -DE-0.3, '$S={}, L={}, J={}, g={}$'.format(s2, l2, j2, np.round(g2, 3)), horizontalalignment='center',verticalalignment='top')    
    
    #plot split d 
    plt.annotate('',((n-0.3)*dx, DE+de1s[-1]), xytext=((n-0.3)*dx, DE+de1s[-2]),
                 arrowprops=dict(facecolor='black',arrowstyle='<->'))
    plt.text((n-0.1)*dx, DE+de1s[-1]-g1/2,'${}\\mu_BB$'.format(np.round(g1, 3)), verticalalignment='center')
    plt.annotate('',((n-0.3)*dx, -DE+de2s[-1]), xytext=((n-0.3)*dx, -DE+de2s[-2]),
                 arrowprops=dict(facecolor='black',arrowstyle='<->'))
    plt.text((n-0.1)*dx, -DE+de2s[-1]-g2/2,'${}\\mu_BB$'.format(np.round(g2, 3)), verticalalignment='center')
    
    plt.axis('off')
    
    #spectrum
    if spectrum:
        low-=1
        d_nu_s, idx = np.unique(d_nu_s, return_index=True)
#        d_nu_s = d_nu_s[idx]
        polars = np.array(polars)
        polars = polars[idx]
        nline = len(d_nu_s)
        line_l = n/2*dx-nline*dx/2
        line_m = n/2*dx
        plt.hlines(low-4, line_l, line_l+nline*dx)
#        for i, d_nu in zip(range(nline),d_nu_s):
#        linexs = line_l+dx/2+dx*np.arange(nline)
        linexs = line_m+d_nu_s/(np.max(d_nu_s)-np.min(d_nu_s))*(nline-1)*dx
        plt.vlines(linexs, low-4, low-2)
        for x, polar, d_nu in zip(linexs, polars, d_nu_s):
            plt.text(x, low-2+0.3, polar,horizontalalignment='center')
            plt.text(x, low-4-0.3, np.round(d_nu, 3), horizontalalignment='center', verticalalignment='top')
        plt.text(line_l-dx, low-4-0.3, '$\\Delta\\tilde{\\nu}/\\tilde{L}$', horizontalalignment='center', verticalalignment='top')
    
        plt.ylim(low-5, high)
    
#    plt.subplot(212)
#    for de1, m1 in zip(de1s, m1s):
        
    
def prompt():
#    detail = False
    detail = input('detail?>>>')
    detail = True if detail in ['1','yes','y','True','true'] else False
#    default = input('default: )
#    if detail in ['1','yes','y','True','true']:
#        s1, l1, j1
    
    s1 = input('Input S1 or "exit" or "set".\nS1 = ')
#    if s1 == 'set':
#        print("detail={}".format(detail))        
    l1 = input('L1 = ')
    j1 = input('J1 = ')
    s2 = input('S2 = ')
    l2 = input('L2 = ')
    j2 = input('J2 = ')
    try: 
        plot(s1, l1, j1, s2, l2, j2)
    except Exception as error:
        print('error')

if __name__ == '__main__':
#    prompt()
#    plot(1,0,1,1,1,2, name1=r'$\mathrm{6s7s{}^3S_1}$', name2=r'$\mathrm{6s6p{}^3P_2}$')
#    plt.figure()
#    plot(0,2,2,0,1,1, name1=r'$\mathrm{6s6d{}^1D_2}$', name2=r'$\mathrm{6s6p{}^1P_1}$')
#    plot(1/2, 1, 1/2, 1/2, 0, 1/2, name1=r'${}^2P_{1/2}$', name2=r'${}^2S_{1/2}$')
    plot(1/2, 1, 3/2, 1/2, 0, 1/2, name1=r'${}^2P_{3/2}$', name2=r'${}^2S_{1/2}$')
    