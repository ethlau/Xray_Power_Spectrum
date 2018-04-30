import numpy as np
import numpy.ma as ma
from scipy import fftpack
from astropy.io import fits

ecf = 1.47e11

n = 17.0
# degree to radian
theta = np.linspace(0., 17/(n-1.0)*360.0)
#
arcsec2rad = np.pi/180.0/3600.0
pixel = 1.2 # in arcsec
tpixel = pixel*arcsec2rad

x = 1.0*np.sin(theta)
y = 1.0*np.cos(theta)

# Fourier space binning
tet_1grid = np.append((1.2+10**(0.15*np.arange(21.0))), 1200.0)
nq_1grid = tet_1grid.size
q_p_minmax = 2.0*np.pi/tet_1grid
k_1_minmax = q_p_minmax/2.0/np.pi
q_p_1grid = np.zeros(nq_1grid-1)
for ik in np.arange(nq_1grid-1) :
    q_p_1grid[ik] = 0.5 * (q_p_minmax[ik]+q_p_minmax[ik+1])

def read_fits_file ( filename, show_hdu_info=False ):

    data, hdr = fits.getdata(filename, header = True)
    return hdr, data;

def dist1D(nrows):
    """Returns a 1-D array in which the value of each element is proportional to its frequency.
    """
    #result = np.linspace(-nrows/2+1, nrows/2, nrows)
    result = fftpack.fftfreq(nrows, d = 1.0/float(nrows))
    #result = fftpack.fftshift(result)
    result = np.sqrt(result**2)
    return result

def calc_del_flux (flux, weight) :

    del_flux = flux - flux.mean()
    del_flux *= weight
    del_flux = del_flux - del_flux.mean()

    return del_flux

def calc_k_w (nx_tot, ny_tot) :
    k_w = np.zeros([nx_tot,ny_tot])
    k_x = dist1D(nx_tot)
    k_y = dist1D(ny_tot)

    for i in np.arange(nx_tot) :
        k_w[i,:] = np.sqrt((k_x[i]/nx_tot)**2.0+(k_y[:]/ny_tot)**2.0)/pixel

    return k_w

def calc_f_clip (cmask, nx_tot, ny_tot) :

    hn0 = np.where(cmask != 0)
    n1 = len(hn0)

    f_clip = float(n1)/nx_tot/ny_tot

    return f_clip

def auto_power ( flux, fab, t_ch, cmask, f_clip, k_w, outfile = None, writefits = None ) :

    s_e = flux.shape
    nx_tot = s_e[0]
    ny_tot = s_e[1]
    nx21 = nx_tot/2+1
    ny21 = ny_tot/2+1

    area = float(nx_tot)*float(ny_tot)*(pixel/3600.0)**2.0/3282.8

    mask = (cmask != 0)

    weight = np.zeros([nx_tot,ny_tot])
    weight = cmask
    #weight_masked = ma.masked_array(weight, ~mask, filled_value = 0.0)

    #print weight_masked.shape, t_ch.shape
    #t_ch_masked = ma.masked_array(t_ch, ~mask)
    #weight_masked = t_ch_masked / t_ch_masked.mean()

    del_flux = calc_del_flux (flux, weight)
    del_flux_ab = calc_del_flux (fab, weight)

    # Remove axis in Fourier Space and compute fft

    #msk_fft = np.ones([nx_tot,ny_tot])
    #msk_fft[:,0] = 0.0
    #msk_fft[0,:] = 0.0
    #msk_fft = np.roll(msk_fft,-nx21, axis = 0)
    #msk_fft = np.roll(msk_fft,-ny21, axis = 1)

    del_flux_fft = fftpack.fft2(del_flux)
    del_flux_ab_fft = fftpack.fft2(del_flux_ab)

    amp = del_flux_fft
    ampab = del_flux_ab_fft
    #amp = fftpack.fftshift(del_flux_fft)
    #ampab= fftpack.fftshift(del_flux_ab_fft)
    #amp = np.roll(np.roll(del_flux_fft, -ny21, axis=1), -nx21, axis = 0)
    #ampab = np.roll(np.roll(del_flux_ab_fft, -ny21, axis=1), -nx21, axis = 0)

    ###
    pairs = np.zeros(nq_1grid-1)
    power = np.zeros(nq_1grid-1)
    sig_p = np.zeros(nq_1grid-1)

    powerab = np.zeros(nq_1grid-1)
    sig_ab = np.zeros(nq_1grid-1)

    k_1_minmax = q_p_minmax/2.0/np.pi
 
    for iq in np.arange(nq_1grid-1) :
        hp = np.where(((k_w >= k_1_minmax[iq+1]) & (k_w < k_1_minmax[iq])))
        if (len(hp) >1) :
            power[iq] = np.mean(np.abs(amp[hp])**2.0)*area/f_clip
            pairs[iq] = len(hp)
            powerab[iq] = np.mean(np.abs(ampab[hp])**2.0)*area/f_clip

    sig_p = power/(np.sqrt(0.5*pairs))
    sig_pab = powerab/(np.sqrt(0.5*pairs))

    # computing P(A+B) - P(A-B)
    pclean = power - powerab
    sig_plc = np.sqrt(sig_pab**2.0 + sig_p**2.0)

    if outfile :
        f = open(outfile,'w')
        for iq in np.arange(nq_1grid-1) :
            print>>f,  2.*np.pi/q_p_1grid[iq], pairs[iq]

        f.close()

    return (amp, ampab, power, powerab, pairs, pclean, sig_plc )

def cross_power ( amp1, amp2, power1, power2, f_clip, k_w, writefits = None ) :

    pairs = np.zeros(nq_1grid-1)
    power = np.zeros(nq_1grid-1)
    sig_p = np.zeros(nq_1grid-1)

    amp_x = amp1.real* amp2.real + amp1.imag * amp2.imag

    for iq in np.arange(nq_1grid -1) :
        hp = where (k_w >= k_1_minmax[iq+1] & k_w < k_1_minmax[iq])

        if (len(hp) >1) :
            power_x[iq] = mean(amp_x[hp])*area/f_clip
            pairs_x[iq] = len(hp)

    sig_x = np.sqrt (0.5*power1*power2/(0.5*pairs_x))

    if writefits :
        pyfits.writeto(writefits, amp_x*area/f_clip)

    return (power_x, sig_x)

def coherence (power_x,pclean1,pclean2,sig_pcl1,sig_pcl2) :

    '''
    !p.multi=[0]
    !p.charsize=1.
    !p.symsize=.1
    SET_PLOT, 'PS'
    DEVICE, FILE='C1X.ps', /COLOR, BITS=8,SCALE_FACTOR=2, YSIZE=13, Xsize=13,/landsc
    ;window, 4, xsize=500,ysize=400

    cirx1=(power1x*power1x)/(pcleanir1*pclean)
    a=(2*power1x/(pcleanir1*pclean)*sig_m1x)^2
    b=((power1x^2/pcleanir1*alog(abs(pclean)))*sig_pcl)^2
    c=((power1x^2/pclean*alog(abs(pcleanir1)))*sig_plcir1)^2
    pippo=sqrt(a+b+c)
    sig_c1=pippo;cirx1*sqrt(12/pairs1x)
    plot,2*!pi/q_p_1grid,cirx1,psym=sym(1),xran=[10,500],yrange=[.001,10],/xstyle,/xlog,/ylog,xtitle=TEXTOIDL('2\pi\q'),ytitle=TEXTOIDL('C_{1,X_{1}}(q)')
    xyouts, 20, 0.800, 'Ch1 vs 0.5-2 keV'
    errplot,2*!pi/q_p_1grid,cirx1-sig_c1,cirx1+sig_c1
    DEVICE, /close
    '''

    return 0

def plot_ps ( power, k, plotfile=None ) :

    import matplotlib.pyplot as plt
    plt.loglog(k, power)
    plt.xlabel(r'$k$')
    plt.ylabel(r'$P(k)$')
    
    if plotfile :
        plt.savefig(plotfile)
    else:
        plt.show()

def main():
    folder = '../Chandra/'
    flux_hdr,flux = read_fits_file (folder+'deltaf.fits')

    t_x_hdr, t_x = read_fits_file (folder+'exp.fits')

    xab_hdr, xab = read_fits_file(folder+'deltaf_ab.fits')

    print ("The mean flux is %e " % ( np.mean(flux) ))

    #cmask = t_ch1 / t_ch1 #;(readfits(folder+''mask_rep_corr.fits',hea5));[0:1131,*]
    cmask_hdr, cmask = read_fits_file(folder+'mask.fits')
    cmask = ( cmask == 1.0 ).astype(int)

    # Fixing the mask after projection
    flux *= cmask/ecf
    xab *= cmask/ecf
    print ("The mean flux after masking is %e " % ( np.mean(flux) ))
    print ("flux array size =", flux.shape)

    s_e = flux.shape
    nx_tot = s_e[0]
    ny_tot = s_e[1]
    nx21 = nx_tot/2+1
    ny21 = ny_tot/2+1

    flux = flux/(pixel/3600.0)**2.0*3282.8 # from cts/s/pix to cts/s/steradians
    xab = xab/(pixel/3600.0)**2.0*3282.8

    # Defining binning in Real and Fourier space
    #exp = np.zeros([nx_tot,ny_tot])
    #rx = np.zeros([nx_tot,ny_tot])
    #n_k = 1000
    #nx_center = nx21 -1
    #ny_center = ny21 -1

    k_w = calc_k_w (nx_tot, ny_tot)
    k_0 = 1.0 / (np.sqrt(2.0) * nx_tot * pixel)
    k_f = 1.01 * np.max(k_w)

    f_clip = calc_f_clip ( cmask, nx_tot, ny_tot )

    (amp, ampab, power, power_ab, pairs, pclean, sig_plc ) = auto_power ( flux, xab, t_x, cmask, f_clip, k_w, outfile = None, writefits = None )

    k_p_1grid = 2.*np.pi/q_p_1grid

    plot_ps(k_p_1grid, pclean)

if __name__ == "__main__":
    main()
