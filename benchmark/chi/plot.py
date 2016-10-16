from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from itertools import product

from params import *

arch = HDFArchive(chi_filename,'r')
ed_arch = HDFArchive(chi_ed_filename,'r')

pp = PdfPages('chi.pdf')
spin_names = ('up','dn')
spin_labels = {'up':'\uparrow\uparrow', 'dn':'\downarrow\downarrow'}

# G_tau
delta_max_text = ""
for sn in spin_names:
    g = arch['G_tau'][sn]
    g_ed = ed_arch['G_tau'][sn]
    oplot(g,    mode='R', lw=0.5, label="QMC, $%s$" % spin_labels[sn])
    oplot(g_ed, mode='R', lw=0.5, label="ED, $%s$" % spin_labels[sn])
    delta_max = np.max(np.abs(g.data[:,0,0] - g_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%s} = %f$\n" % (spin_labels[sn],delta_max)
ax = plt.gca()
ax.set_title('$G(\\tau)$')
ax.set_ylabel('$G(\\tau)$')
ax.set_xlim((0,beta))
ax.set_ylim((-1,0.05))
ax.legend(loc='lower center',prop={'size':10})
ax.text(beta/2,-0.5,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# G_iw
plt.cla()
delta_max_text = ""
for sn in spin_names:
    g = arch['G_iw'][sn]
    g_ed = ed_arch['G_iw'][sn]
    oplot(g,    mode='I', lw=0.5, label="QMC, $%s$" % spin_labels[sn])
    oplot(g_ed, mode='I', lw=0.5, label="ED, $%s$" % spin_labels[sn])
    delta_max = np.max(np.abs(g.data[:,0,0] - g_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%s} = %f$\n" % (spin_labels[sn],delta_max)
ax = plt.gca()
ax.set_title('$G(i\\omega)$')
ax.set_ylabel('$G(i\\omega)$')
ax.set_xlim((0,5.0))
ax.legend(loc='upper center',prop={'size':10})
ax.text(2.5,0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# G_l
plt.cla()
delta_max_text = ""
for sn in spin_names:
    g = arch['G_l'][sn]
    g_ed = ed_arch['G_l'][sn]
    oplot(g,    mode='R', lw=0.5, label="QMC, $%s$" % spin_labels[sn])
    oplot(g_ed, mode='R', lw=0.5, label="ED, $%s$" % spin_labels[sn])
    delta_max = np.max(np.abs(g.data[:,0,0] - g_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%s} = %f$\n" % (spin_labels[sn],delta_max)
ax = plt.gca()
ax.set_title('$G(\\ell)$')
ax.set_ylabel('$G(\\ell)$')
ax.set_xlim((0,n_l-1))
ax.legend(loc='lower center',prop={'size':10})
ax.text(n_l/2,-2,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# chi_tau (diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    chi = arch['chi_tau'][s,s]
    chi_ed = ed_arch['chi_tau'][s,s]
    oplot(chi,    mode='R', lw=0.5, label="QMC, %i%i" % (s,s))
    oplot(chi_ed, mode='R', lw=0.5, label="ED, %i%i" % (s,s))
    delta_max = np.max(np.abs(chi.data[:,0,0] - chi_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(\\tau)$ (diagonal)')
ax.set_ylabel('$\\chi(\\tau)$')
ax.set_xlim((0,beta))
ax.set_ylim((0,1.0))
ax.legend(loc='upper center',prop={'size':10})
ax.text(beta/2,0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# chi_tau (off-diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    chi = arch['chi_tau'][s,1-s]
    chi_ed = ed_arch['chi_tau'][s,1-s]
    oplot(chi,    mode='R', lw=0.5, label="QMC, %i%i" % (s,1-s))
    oplot(chi_ed, mode='R', lw=0.5, label="ED, %i%i" % (s,1-s))
    delta_max = np.max(np.abs(chi.data[:,0,0] - chi_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,1-s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(\\tau)$ (off-diagonal)')
ax.set_ylabel('$\\chi(\\tau)$')
ax.set_xlim((0,beta))
ax.set_ylim((0,0.5))
ax.legend(loc='lower center',prop={'size':10})
ax.text(beta/2,0.2,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# chi_iw (diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    chi = arch['chi_iw'][s,s]
    chi_ed = ed_arch['chi_iw'][s,s]
    oplot(chi,    mode='R', lw=0.5, label="QMC, %i%i" % (s,s))
    oplot(chi_ed, mode='R', lw=0.5, label="ED, %i%i" % (s,s))
    delta_max = np.max(np.abs(chi.data[:,0,0] - chi_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(i\\omega)$ (diagonal)')
ax.set_ylabel('$\\chi(i\\omega)$')
ax.set_xlim((0,5.0))
ax.legend(loc='upper center',prop={'size':10})
ax.text(2.5,ax.get_ylim()[1]*0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# chi_iw (off-diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    chi = arch['chi_iw'][s,1-s]
    chi_ed = ed_arch['chi_iw'][s,1-s]
    oplot(chi,    mode='R', lw=0.5, label="QMC, %i%i" % (s,1-s))
    oplot(chi_ed, mode='R', lw=0.5, label="ED, %i%i" % (s,1-s))
    delta_max = np.max(np.abs(chi.data[:,0,0] - chi_ed.data[:,0,0]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,1-s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(i\\omega)$ (off-diagonal)')
ax.set_ylabel('$\\chi(i\\omega)$')
ax.set_xlim((0,5.0))
ax.legend(loc='upper center',prop={'size':10})
ax.text(2.5,ax.get_ylim()[1]*0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

chi_l_ed = ed_arch['chi_l']
chi_l = chi_l_ed.copy()
chi_l << MatsubaraToLegendre(arch['chi_iw'])

# chi_l (diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    oplot(chi_l[s,s],    mode='R', lw=0.5, label="QMC, $%i%i$" % (s,s))
    oplot(chi_l_ed[s,s], mode='R', lw=0.5, label="ED, $%i%i$" % (s,s))
    delta_max = np.max(np.abs(chi_l.data[:,s,s] - chi_l_ed.data[:,s,s]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(\\ell)$ (diagonal)')
ax.set_ylabel('$\\chi(\\ell)$')
ax.set_xlim((0,n_l-1))
ax.legend(loc='upper center',prop={'size':10})
ax.text(n_l/2,ax.get_ylim()[1]*0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

# chi_l (off-diagonal)
plt.cla()
delta_max_text = ""
for s in (0,1):
    oplot(chi_l[s,1-s],    mode='R', lw=0.5, label="QMC, $%i%i$" % (s,1-s))
    oplot(chi_l_ed[s,1-s], mode='R', lw=0.5, label="ED, $%i%i$" % (s,1-s))
    delta_max = np.max(np.abs(chi_l.data[:,s,1-s] - chi_l_ed.data[:,s,1-s]))
    delta_max_text += "$\\delta^{max}_{%i%i} = %f$\n" % (s,1-s,delta_max)
ax = plt.gca()
ax.set_title('$\\chi(\\ell)$ (off-diagonal)')
ax.set_ylabel('$\\chi(\\ell)$')
ax.set_xlim((0,n_l-1))
ax.legend(loc='upper center',prop={'size':10})
ax.text(n_l/2,ax.get_ylim()[1]*0.6,delta_max_text,horizontalalignment='center')
pp.savefig(plt.gcf())

pp.close()
