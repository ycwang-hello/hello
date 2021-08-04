# -*- coding: utf-8 -*-
"""
Created on Sun Aug  1 15:07:30 2021

@author: Yuchen Wang

可以修改217行前后的beta变量观察光源速度不同时的情况
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import warnings

def sr_v_inv_trans2(u):
    '''
    2-d special relativity inverse velocity transforamtion
    u: Sigma'系在Sigma系中的，x方向牵连速度/光速 （即几何单位制下的速度，或者beta）
    在Sigma'系中发光，要变换到Sigma系，的洛伦兹变换。
    '''
    assert u <= 1, '超光速不允许！'
    gamma = 1/np.sqrt(1-u**2)
    def transform(vxp, vyp):
        '''
        vxp, vyp: 在Sigma'系中的速度/光速，或者说几何单位制下的数。
        '''
        assert np.max(np.sqrt(vxp**2+vyp**2)) <= 1, '超光速不允许！'
        vy = vyp/(gamma*(1+u*vxp))
        vx = (vxp+u)/(1+u*vxp)
        return vx, vy
    return transform

class WaveFront():

    # def ellipse(a, x0, y0, e=0):
    #     thetas = np.linspace(0, 2*np.pi, 300)
    #     c = a * e#np.sqrt(a**2-b**2)
    #     b = np.sqrt(a**2 - c**2)
    #     x = a*np.cos(thetas) + x0
    #     y = b*np.sin(thetas) + y0
    #     return x, y, c
    def __init__(self, x, y, c, v_src=None, t0=0, life=np.inf, wf_n=15, color='black', show_v=True, ax=None):
        '''
        
        Parameters
        ----------
        x : initial coordinate x.
        y : initial coordinate x.
        c : velocity of the wave.
        wf_n : number of points on the wavefront.
        life: the life span of a wavefront.
        color : color of the points on the wavefront.
        ax : matplotlib.axes._subplots.AxesSubplot
            axis

        '''
        self.x = x
        self.y = y
        self.c = c
        self.ax = ax
        r0 = c*t0
        self.life = life
        self.t = t0
        self.dead = False
        if v_src == None:
            warnings.warn('未输入波源速度，无法考虑相对论集束效应。')
            v_src = 0
        self.v_inv_tran = sr_v_inv_trans2(v_src)
        if self.ax == None:
            self.ax = plt.gca()
        self.r = r0
        self.wf_n = wf_n
        self.color = color
        self.show_v = show_v
        # xs, ys, _ = WaveFront.ellipse(self.r, self.x, self.y)
        
        # plot points on the wavefront
        self.xs, self.ys = np.ones(self.wf_n)*self.x, np.ones(self.wf_n)*self.y
        self.thetas = np.linspace(0, 2*np.pi, self.wf_n)
        self.vxps, self.vyps = self.c * np.cos(self.thetas), self.c * np.sin(self.thetas)
        self.vxs, self.vys = self.v_inv_tran(self.vxps, self.vyps)
        self.wf_pt, = self.ax.plot(self.xs, self.ys, linestyle='', marker='.', markersize=7, color=self.color)
        
        # plot arrows showing the velocity of points (photons) on the wavefront
        if self.show_v:
            self.wf_arrs = []
            for xsi, ysi, vxsi, vysi in zip(self.xs, self.ys, self.vxs, self.vys):
                arr = self.ax.arrow(xsi, ysi, vxsi, vysi)
                self.wf_arrs.append(arr)
        
        # plot the wavefront
        self.lxs, self.lys = np.ones(300)*self.x, np.ones(300)*self.y
        self.lthetas = np.linspace(0, 2*np.pi, 300)
        self.lvxps, self.lvyps = self.c * np.cos(self.lthetas), self.c * np.sin(self.lthetas)
        self.lvxs, self.lvys = self.v_inv_tran(self.lvxps, self.lvyps)
        self.wf_line, = self.ax.plot(self.lxs, self.lys, 'k--')
        
        self.update(t0)
        
    def update(self, dt=.1):
        if self.t >= self.life and not self.dead:
            self.dead = True
            # self.circle.remove()
            # self.circle.set_linestyle('')
            self.wf_line.remove()
            self.wf_pt.remove()
            if self.show_v:
                for arr in self.wf_arrs:
                    arr.remove()
        if not self.dead:
            self.r += self.c * dt 
            self.t += dt
            # xs, ys, _ = WaveFront.ellipse(self.r, self.x, self.y)
            self.xs += self.vxs*dt
            self.ys += self.vys*dt
            self.wf_pt.set_xdata(self.xs)
            self.wf_pt.set_ydata(self.ys)

            self.lxs += self.lvxs*dt
            self.lys += self.lvys*dt
            self.wf_line.set_xdata(self.lxs)
            self.wf_line.set_ydata(self.lys)
            
            if self.show_v:
                for arr in self.wf_arrs:
                    arr.remove()
            
                self.wf_arrs = []
                for xsi, ysi, vxsi, vysi in zip(self.xs, self.ys, self.vxs, self.vys):
                    arr = self.ax.arrow(xsi, ysi, vxsi, vysi)
                    self.wf_arrs.append(arr)


class Source():
    txt_disp = .2 #text displacement
    def __init__(self, x, y, v, c, f, theta=0, wflife=np.inf, ax=None, autoview=False, show_wf_v=True):
        self.x = x
        self.y = y
        self.v = v # the velocity of source
        self.theta = theta # the direction of velocity of source
        self.vx = v*np.cos(theta)
        self.vy = v*np.sin(theta)
        self.c = c # the velocity of wave
        self.f = f
        self.u = lambda t: np.cos(2*np.pi*self.f*t) #u is the displacement
        self.ax = ax
        self.t = 0
        self.wflife = wflife #life of wavefronts
        self.p = 1/f # p is period
        self.autoview = autoview
        self.show_wf_v = show_wf_v
        if self.ax == None:
            self.ax = plt.gca()
            
        self.src_pt, = plt.plot([x], [y], marker='*', markersize=10)
        self.src_txt = plt.text(x, y-Source.txt_disp, 'src', fontsize=12, horizontalalignment='center', verticalalignment='top')
            
        self.wavefronts = []
        self.wavefronts.append(WaveFront(self.x, self.y, self.c, v_src=self.v, life=self.wflife, ax=self.ax, show_v=self.show_wf_v))
        
    def update(self, dt=.1):
        assert dt < self.p
        
        # wavefront propagate
        for i, wavefront in enumerate(self.wavefronts):
            wavefront.update(dt)
                
        # source move
        self.x += self.vx * dt
        self.y += self.vy * dt
        self.src_pt.set_xdata([self.x])
        self.src_pt.set_ydata([self.y])
        self.src_txt.set_x(self.x)
        self.src_txt.set_y(self.y-Source.txt_disp)
        
        # generate new wavefront if needed
        if (self.t+dt) // self.p - self.t // self.p == 1:
            t0 = self.t+dt - ((self.t+dt)//self.p)*self.p
            x0 = self.x - self.vx * t0
            y0 = self.y - self.vy * t0
            self.wavefronts.append(WaveFront(x0, y0, self.c, t0=t0, v_src=self.v, life=self.wflife, ax=self.ax, show_v=self.show_wf_v))
        self.t += dt    
        
        if self.autoview:
            # recompute the ax.dataLim
            self.ax.relim()
            # update ax.viewLim using the new dataLim
            self.ax.autoscale_view()

# dt = .05
# ts = np.arange(0, 9.5, dt)
# fig = plt.figure(figsize=(6, 6))
# ax = fig.subplots()
# ax.set_aspect('equal')
# ax.set_xlim([-10, 13])
# ax.set_ylim([-10, 10])
# src = Source(x=0, y=0, v=.9, c=1, f=.3, wflife=12, ax=ax, theta=0)#np.pi/4)
# def update(t):
#     src.update(dt)
#     src.ax.set_title('t={:.2f}'.format(t))

# ani = animation.FuncAnimation(
#     fig, update, frames=ts, interval=.5*dt*1e3, 
#     blit=False, save_count=len(ts), repeat=False) 

# t=0 时，光源应该在原点
show_k = True
show_kp = True
show_phi = False
show_phase = False #False  #等相位面
show_realdir = True
show_z = True
show_z_cont_label = False #show z contour label on the lines

isophase_n = 10

beta = .7 #.9#.999#.5#.9
gamma = 1/np.sqrt(1-beta**2)
freq = .5 #实际是freq/c
c = 1

dct = .05
cts = np.arange(0, 15, dct)
fig = plt.figure(figsize=(10, 8))
ax = fig.subplots()
ax.set_aspect('equal')
ax.set_xlim([-10, 13])
ax.set_ylim([-10, 10])

x = np.linspace(-13, 13, 100)
y = x
app = 10 
X, Y = np.meshgrid(x, y)

dx = (x[1]-x[0])/2.
dy = (y[1]-y[0])/2.
extent = [x[0]-dx, x[-1]+dx, y[0]-dy, y[-1]+dy]

cos_phi = X/np.sqrt(X**2+Y**2)
sin_phi = Y/np.sqrt(X**2+Y**2)

cos_realang = lambda ct: (X-beta*ct)/np.sqrt((X-beta*ct)**2+Y**2)
sin_realang = lambda ct: Y/np.sqrt((X-beta*ct)**2+Y**2)

cos_theta_p = lambda ct: gamma*(X-beta*ct)/np.sqrt(gamma**2*(X-beta*ct)**2+Y**2)
sin_theta_p = lambda ct: Y/np.sqrt(gamma**2*(X-beta*ct)**2+Y**2)
omega_over_omegap = lambda ct: gamma*(1+beta*cos_theta_p(ct))
cos_theta = lambda ct: 1/omega_over_omegap(ct)*gamma*(beta+cos_theta_p(ct))
sin_theta = lambda ct: 1/omega_over_omegap(ct)*sin_theta_p(ct)

phase = lambda ct: omega_over_omegap(ct)*(ct - X*cos_theta(ct) - Y*sin_theta(ct))
z = lambda ct: 1/omega_over_omegap(ct) - 1
z = lambda ct: -np.log(omega_over_omegap(ct)) #以下名为z的实际上是ln(1+z)
phi0 = np.arctan2(np.sqrt(2*gamma**2), -np.sqrt(gamma-1))

zmin, zmax = np.min(z(0)), np.max(z(0))
n_colorbar_ticks = 8
ticks = np.arange(n_colorbar_ticks+1)
ticks = ticks / n_colorbar_ticks * (zmax-zmin) + zmin #此处ticks为z
ticklabels = np.exp(ticks) - 1
ticklabels += 5e-5
ticklabels = ['{:.2f}'.format(ticklabel) for ticklabel in ticklabels]
ticks *= -1 #画的是-z，所以tick取负号
lines = ticks.copy()
line_per_clabel = 3

def get_z0line(ct):
    l = 20
    x0 = beta*ct
    xs = [x0+l*np.cos(phi0), x0, x0+l*np.cos(phi0)]
    ys = [0+l*np.sin(phi0), 0, 0-l*np.sin(phi0)]
    return xs, ys

contourf_kwargs = dict(alpha=.5, cmap='RdBu')
contour_kwargs = dict(alpha=.8, cmap='RdBu', levels=20, linestyles='dashed') #, levels=20, label='$\\ln(1+z)$'
heatmap_kwargs = dict(alpha=.5, cmap='RdBu', extent=extent, interpolation='bilinear')#, label='$\\ln(1+z)$')
z_clabel_kwargs = dict(fontsize=10, colors='k', fmt='%1.2f')

if show_z:
    z_contour = ax.contour(X, Y, -z(0), **contour_kwargs)
    z_heatmap = ax.imshow(-z(0), **heatmap_kwargs)
    xs, ys = get_z0line(0)
    z_eq0, = ax.plot(xs, ys, 'k-.', label='$z=0$')
    z_colorbar = fig.colorbar(z_heatmap, ax=ax, shrink=.7, label='$z$', ticks=ticks) # label='$\\ln(1+z)$'
    z_colorbar.ax.set_yticklabels(ticklabels)
    if show_z_cont_label:
        z_clabel = z_contour.clabel(z_contour.levels[::line_per_clabel], **z_clabel_kwargs) #add labels to the contour
if show_kp:
    kp = ax.quiver(X[::app, ::app], Y[::app, ::app], 
                   cos_theta_p(0)[::app, ::app], sin_theta_p(0)[::app, ::app],
                   color='blue', label='$\\vec{k}\'$')#, scale=.01)
if show_phase:
    isophase = ax.contour(X, Y, phase(0), levels=isophase_n, label='isophase')
if show_phi:
    phi_arrow = ax.quiver(X[::app, ::app], Y[::app, ::app], 
                          cos_phi[::app, ::app], sin_phi[::app, ::app],
                          color='black', label='$\\hat{r}$')#, scale=.01)
if show_realdir:
    real_ang = ax.quiver(X[::app, ::app], Y[::app, ::app], 
                         cos_realang(0)[::app, ::app], sin_realang(0)[::app, ::app],
                         color='black', label='$\\vec{r} - \\vec{r}_{\\mathrm{src}}$')#, scale=.01)
if show_k:
    k = ax.quiver(X[::app, ::app], Y[::app, ::app], 
                  cos_theta(0)[::app, ::app], sin_theta(0)[::app, ::app],
                  color='red', label='$\\vec{k}$')#, scale=.01)

ax.legend()
src = Source(x=0, y=0, v=beta*c, c=c, f=freq*c, wflife=12, ax=ax, theta=0, show_wf_v=False)#np.pi/4)

def update(ct):
    vt = beta*ct
    global isophase, z_contour, z_clabel
    src.update(dct)
    src.ax.set_title('$\\beta={:.2f},\\ ct={:.2f}$'.format(beta, ct))
    if show_z:
        for coll in z_contour.collections:
            ax.collections.remove(coll)

        z_contour = ax.contour(X, Y, -z(ct), **contour_kwargs)
        z_heatmap.set_data(-z(ct))
        xs, ys = get_z0line(ct)
        z_eq0.set(xdata=xs, ydata=ys)
        if show_z_cont_label:
            for i in range(len(z_clabel)):
                z_clabel[i].remove()
            z_clabel = z_contour.clabel(z_contour.levels[::line_per_clabel], **z_clabel_kwargs) #add labels to the contour

    if show_kp:
        kp.set_UVC(cos_theta_p(ct)[::app, ::app], sin_theta_p(ct)[::app, ::app])
    if show_k:
        k.set_UVC(cos_theta(ct)[::app, ::app], sin_theta(ct)[::app, ::app])
    if show_phase:
        for coll in isophase.collections:
            ax.collections.remove(coll)
        isophase = ax.contour(X, Y, phase(ct), levels=isophase_n, label='isophase')

    if show_phi:
        pass #phi_arrow.set_UVC(cos_phi[::app, ::app], sin_phi[::app, ::app])
    if show_realdir:
        real_ang.set_UVC(cos_realang(ct)[::app, ::app], sin_realang(ct)[::app, ::app])

plt.tight_layout()

ani = animation.FuncAnimation(
    fig, update, frames=cts, interval=dct*1e3, 
    blit=False, save_count=5,#len(cts), 
    repeat=False,
    ) 

ani.save("SR_Dopller_aberration_full_1.mp4")

