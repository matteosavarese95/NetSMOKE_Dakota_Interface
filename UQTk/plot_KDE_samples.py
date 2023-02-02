import sys
import os

sys.path.append(os.environ["UQTK_INS"]+'/PyUQTk')

from utils import pdf_kde
    
try:
    from scipy import stats
except ImportError:
    print("Scipy stats module could not be found")
    
from datetime import datetime    
import matplotlib.pyplot as plt   
import numpy as np
plt.style.use('bello.mplstyle') #reads matplotlib custom style
plt.rcParams['text.usetex'] = True
        
###################################################################################
def KDE(fcn_evals):
    """
    Performs kernel density estimation
    Input:
        fcn_evals: numpy array of evaluations of the forward model (values of heat flux Q)
    Output:
        xpts_pce: numpy array of points at which the PDF is estimated.
        PDF_data_pce: numpy array of estimated PDF values.
    """
    # Perform KDE on fcn_evals
    kern_pce=stats.kde.gaussian_kde(fcn_evals)
    # Generate points at which to evaluate the PDF
    xpts=np.linspace(fcn_evals.min(),fcn_evals.max(),200)
    # Evaluate the estimated PDF at these points
    PDF_data=kern_pce(xpts)
    return xpts, PDF_data
    

def plot_pdf2d(samsx, samsy, pltype='kde', ncont=10, color=None, lwidth=1, mstyle='o', ax=None, limits=None):
    if ax is None:
        ax = plt.gca()

    if pltype == 'kde':
        ngrid = 400
        a = 4.
        pnomx = np.mean(samsx)
        pdelx = np.std(samsx)
        pnomy = np.mean(samsy)
        pdely = np.std(samsy)

        x = np.linspace(pnomx - a * pdelx, pnomx + a * pdelx, ngrid)
        y = np.linspace(pnomy - a * pdely, pnomy + a * pdely, ngrid)
        if limits:
            x = np.linspace(limits[0][0], limits[0][1], ngrid)
            y = np.linspace(limits[1][0], limits[1][1], ngrid)
        X, Y = np.meshgrid(x, y)
        pgrid = np.vstack((X.flatten(), Y.flatten())).T  # pgrid.shape is (33^2,2)
        _, pdf = pdf_kde.get_pdf(np.vstack((samsx, samsy)).T, pgrid, verbose=0, method='Python')

        if color is None:
            ax.contour(X, Y, pdf.reshape(X.shape), ncont, linewidths=lwidth)
        else:
            ax.contour(X, Y, pdf.reshape(X.shape), ncont, colors=color, linewidths=lwidth)

    elif pltype == 'sam':
        ax.plot(samsx, samsy, color=color, marker=mstyle, linestyle='None')
    else:
        print("Plot type is not recognized. Exiting")
        sys.exit()

    ax.grid(False)

    return ax


def plot_grid_pdf2d(var1,var2):
    varx = obstitle.index(var1)
    vary = obstitle.index(var2)
    fig,ax = plt.subplots(nphi,nTin,figsize=(10,10))
    fig.subplots_adjust(bottom=0.2)
    for i in range(nphi):
        for j in range(nTin):
            ax[i,j].set_xlim(xlimits[varx])
            ax[i,j].set_ylim(xlimits[vary])
            plot_pdf2d(response[:,nTin*nobs*i+nobs*j+varx],response[:,nTin*nobs*i+nobs*j+vary], ncont = 10, ax=ax[i,j], limits=(xlimits[varx],xlimits[vary]))
            if(i == 0): ax[i,j].set_title(f'T = {xcond[nTin*i+j,1]}')
            if(j == 0): ax[i,j].set_ylabel(r'$\Phi$'+f' = {xcond[nTin*i+j,0]}',rotation=0, labelpad=20)
    fig.supxlabel(obslabel[varx])
    fig.supylabel(obslabel[vary])
    plt.savefig(f"joint_{obstitle[varx]}_{obstitle[vary]}.png",dpi=300, transparent=False)
    plt.close()
    

def plot_grid_pdf2d_corr(var,par):
    varx = obstitle.index(var)
    vary = parstitle.index(par)
    fig,ax = plt.subplots(nphi,nTin,figsize=(10,10))
    fig.subplots_adjust(bottom=0.2)
    for i in range(nphi):
        for j in range(nTin):
            ax[i,j].set_xlim(xlimits[varx])
            ax[i,j].set_ylim(plimits[vary])
            plot_pdf2d(response[:,nTin*nobs*i+nobs*j+varx],inpars[:,vary], ncont = 10, ax=ax[i,j], limits=(xlimits[varx],plimits[vary]))
            if(i == 0): ax[i,j].set_title(f'T = {xcond[nTin*i+j,1]}')
            if(j == 0): ax[i,j].set_ylabel(r'$\Phi$'+f' = {xcond[nTin*i+j,0]}',rotation=0, labelpad=20)
    fig.supxlabel(obslabel[varx])
    fig.supylabel(parslabel[vary])
    plt.savefig(f"joint_{obstitle[varx]}_{parstitle[vary]}.png",dpi=300, transparent=False)
    plt.close()
 
###################################################################################    
now = datetime.now()
date = now.strftime("%Y-%m-%dT%H%M%S")
folder = 'KDE_'+date
os.system(f'mkdir {folder}')

inpars = np.loadtxt('ptrain.dat')
response = np.loadtxt('ytrain.dat')    
xcond = np.loadtxt('xcond.dat')
ncond = xcond.shape[0]
xlabel = ['phi_rich', 'T_in']

npars = 4  #'T_cstr_1', 'TAU_cstr_2','H_cstr_2', 'L_pfr' 
nobs = 4   #'T_cstr_2', 'T_pfr_4', 'NO_pfr_4', 'NH3_pfr_4'


print(f'Found {response.shape[0]} samples')

nphi = np.count_nonzero(xcond == xcond[0,1])
nTin = np.count_nonzero(xcond == xcond[0,0])

print(f'Found {response.shape[1]} responses')
print(f'Found {ncond} x-conditions: {nphi} Phi, {nTin} Tin')

    
obstitle = ['T_cstr_2', 'T_pfr_4', 'NO_pfr_4', 'NH3_pfr_4'] * ncond
obslabel = ['T [K]', 'T [K]', 'NO ppm', 'NH3 ppm'] * ncond
xlimits = [(1900,2300),(1170,1230),(-10,600),(-10,600)] * ncond


parstitle = ['T_cstr_1', 'TAU_cstr_2','H_cstr_2', 'L_pfr'] 
parslabel = ['T [K]', r'\tau_{CSTR2} [s]',r'h_{CSTR2} [J/m^2/K/s]', r'L_{PFR} [mm]'] 
plimits = [(600,950),(0.04,0.06),(0.9411768,1.4117652),(2000,6000)] 

# plot KDE of PDFs of single responses

allxpts = []
allPDF = []
for resp in range(response.shape[1]):
    xpts, PDF_data= KDE(response[:,resp])
    allxpts.append(xpts)
    allPDF.append(PDF_data)


cond = 0   
for resp in range(response.shape[1]):
    if(resp >0 and resp % 4 == 0): cond=cond+1
    #xpts, PDF_data= KDE(response[:,resp])
    fig,ax = plt.subplots(figsize=(4,3))
    plt.plot(allxpts[resp], allPDF[resp], linewidth=2, color='r', label= f'KDE, Monte Carlo sampling')
    plt.title(f'{obstitle[resp]} - phiR {xcond[cond,0]} - Tin {xcond[cond,1]}')
    plt.xlim(xlimits[resp])
    plt.xlabel(f'{obslabel[resp]}')
    plt.savefig(f"{folder}/KDE_{obstitle[resp]}_cond{cond}.png",dpi=800, transparent=False)
    plt.close()
    

# plot joint PDFs

plot_grid_pdf2d('T_cstr_2','NO_pfr_4')

plot_grid_pdf2d('NH3_pfr_4','NO_pfr_4')

plot_grid_pdf2d_corr('NO_pfr_4','T_cstr_1')
