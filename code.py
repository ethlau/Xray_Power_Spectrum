import numpy as np
import numpy.ma as ma
import scipy
import pyfits

ecf = 1.47e11

n = 17.0
# degree to radian
theta = np.linspace(17/(n-1.0)*360.0)
#
arcsec2rad = np.pi/180.0/3600.0
pixel = 1.2 # in arcsec
tpixel = pixel*arcsec2rad

x = 1.0*sin(theta)
y = 1.0*cos(theta)

tet_1grid = np.append((1.2+10**(0.15*np.arange(21.0))), 1200.0)
nq_1grid = tet_1grid.size
q_p_minmax = 2.0*np.pi/tet_1grid
k_1_minmax = q_p_minmax/2.0/np.pi
q_p_1grid = np.zeros(nq_1grid-1)
for ik in np.arange(nq_1grid-1) :
    q_p_1grid[ik] = 0.5 * (q_p_minmax[ik]+q_p_minmax[ik+1])

def read_fits_file ( filename ):

    data = pyfits.getdata( filename )

    return data;

def dist1D(nrows):
    """Returns a 1-D array in which the value of each element is proportional to its frequency.
    """
    result = np.linspace(-nrows/2+1, nrows/2, nrows)
    result = np.sqrt(result**2)
    return result

def calc_del_flux (flux, mask, weight_masked) :

    del_flux_masked = ma.masked_array(flux, ~mask, filled_value = 0.0))
    del_flux_masked = del_flux_masked - del_flux_masked.mean()
    del_flux_masked *= weight_masked
    del_flux_masked = del_flux_masked - del_flux_masked.mean()

    return del_flux_masked

def calc_k_w (nx_tot, ny_tot) :
    k_w = np.zeros([nx_tot,ny_tot])
    k_x = dist1D(nx_tot)
    k_y = dist1D(ny_tot)

    for i in np.arange(nx_tot) :
        k_w[:,i] = np.sqrt((k_x[i]/nx_tot)**2.0+(k_y[:]/ny_tot)**2.0)/pixel

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

    area = float(nx_tot)*float(ny_tot)*/(pixel/3600.0)**2.0/3282.8

    mask = (cmask != 0)

    weight = np.zeros([nx_tot,ny_tot])
    weight_masked = ma.masked_array(weight, ~mask, filled_value = 0.0)
    t_ch_masked = ma.masked_array(t_ch, ~mask)
    weight_masked = t_ch_masked / t_ch_masked.mean()

    del_flux_masked = calc_del_flux (flux, mask, weight_masked)
    del_flux_ab_masked = calc_del_flux (fab, mask, weight_masked)

    # Remove axis in Fourier Space and compute fft

    msk_fft = np.ones([nx_tot,ny_tot])
    msk_fft[:,0] = 0.0
    msk_fft[0,:] = 0.0
    msk_fft = np.roll(msk_fft,-nx21, axis = 0)
    msk_fft = np.roll(msk_fft,-ny21, axis = 1)

    del_flux_fft = np.fft2(del_flux_masked)
    del_flux_ab_fft = np.fft2(del_flux_ab_masked)

    amp = np.roll(np.roll(del_flux_fft, -ny21, axis=1), -nx21, axis = 0)
    ampab = np.roll(np.roll(del_flux_ab_fft, -ny21, axis=1), -nx21, axis = 0)

    ###
    pairs = np.zeros(nq_1grid-1)
    power = np.zeros(nq_1grid-1)
    sig_p = np.zeros(nq_1grid-1)

    powerab = np.zeros(nq_1grid-1)
    sig_ab = np.zeros(nq_1grid-1)

    k_1_minmax = q_p_minmax/2.0/np.pi

    for iq in np.arange(nq_1grid-1) :
        hp = where (k_w >= k_1_minmax[iq+1] & k_w < k_1_minmax[iq])

        if (len(hp) >1) :
            power[iq] = mean(abs(amp[hp])**2.0)*area/f_clip
            pairs[iq] = len(hp)
            powerab[iq] = mean(abs(ampab[hp])**2.0)*area/f_clip

    sig_p = power/(np.sqrt(0.5*pairs))
    sig_pab = powerab/(np.sqrt(0.5*pairs))

    # computing P(A+B) - P(A-B)
    pclean = power - powerab
    sig_plc = np.sqrt(sig_pab**2.0 + sig_p**2.0)

    if outfile :
        for iq in no.arange(nq_1grid-1) :
            print>>outfile,  2.*np.pi/q_p_1grid[iq], pairs[iq]

        outfile.close()

    if writefits :
        pyfits.writeto(writefits,amp.real+amp.imag)

    return ( amp, amp_ab, power, power_ab, pclean, sig_plc )

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

def main():
    ''' Note: difference between IDL and python
        IDL: column major order; python: row major order
    '''

    flux1 = read_fits_file ('egs_flu_ch1.fits')
    fab1  = read_fits_file ('ch1_ab_rep.fits')
    t_ch1 = read_fits_file ('egs_exp_ch1.fits')
    mean2 = read_fits_file ('mean_0520_flux.fits')

    flux2 = read_fits_file ('egs_flu_ch2.fits')
    fab2  = read_fits_file ('ch2_ab_rep.fits')
    t_ch2 = read_fits_file ('egs_exp_ch2.fits')

    flux4 = read_fits_file ('egs_flu_ch4.fits')
    fab4  = read_fits_file ('ch4_ab_rep.fits')
    t_ch4 = read_fits_file ('egs_exp_ch2.fits')

    img = read_fits_file ('flux_rep.fits')

    flux = read_fits_file ('deltaf_0520_new.fits')

    t_x = read_fits_file ('exp_ab_0520.fits')

    cmask = t_ch1 / t_ch1 #;(readfits('mask_rep_corr.fits',hea5));[0:1131,*]
    xab = read_fits_file('deltaf_ab_0520.fits')

    print ("The mean flux is %e " % ( np.mean(flux) ))

    # Fixing the mask after projection
    cmask[cmask < 1.0] = 0.0
    cmask[cmask >= 1.0] = 1.0

    pyfits.writeto('cmask_fixex.futs',cmask)

    flux1 *= mask
    flux2 *= mask
    flux3 *= mask
    flux *= mask/ecf
    mean2 *= mask/ecf
    xab *= mask/ecf

    fab1 *= mask
    fab2 *= mask
    fab4 *= mask

    print ("The mean flux after masking is %e " % ( np.mean(flux) ))

    s_e = flux1.shape
    nx_tot = s_e[0]
    ny_tot = s_e[1]
    nx21 = nx_tot/2+1
    ny21 = ny_tot/2+1

    flux = flux/(pixel/3600.0)**2.0*3282.8 # from cts/s/pix to cts/s/steradians
    xab = xab/(pixel/3600.0)**2.0*3282.8

    # Defining binning in Real and Fourier space
    exp = np.zeros([nx_tot,ny_tot])
    rx = np.zeros([nx_tot,ny_tot])
    n_k = 1000
    nx_center = nx21 -1
    ny_center = ny21 -1

    kw = calc_k_w (nx_tot, ny_tot)
    k_0 = 1.0 / (np.sqrt(2.0) * nx_tot * pixel)
    k_f = 1.01 * np.max(k_w)

    f_clip = calc_f_clip ( cmask, nx_tot, ny_tot )

    # IR Channel 1
    (amp_ir1, amp_ab_ir1, power_ir1, power_ab_ir1, pclean_ir1, sig_plc_ir1 ) = auto_power ( flux1, fab1, t_ch1, cmask, f_clip, k_w, outfile = 'pairs_egs.dat', writefits = None )

    # IR Channel 2
    (amp_ir2, amp_ab_ir2, power_ir2, power_ab_ir2, pclean_ir2, sig_plc_ir2 ) = auto_power ( flux2, fab2, t_ch2, cmask, f_clip, k_w, outfile = None, writefits = None )

    # Channel 12
    (amp_ir4, amp_ab_ir4, power_ir4, power_ab_ir4, pclean_ir4, sig_plc_ir4 ) = auto_power ( flux4, fab4, t_ch4, cmask, f_clip, k_w, outfile = None, writefits = None )

    # X-ray
    (amp, amp_ab, power, power_ab, pclean, sig_plc ) = auto_power ( flux, fab, t_x, cmask, f_clip, k_w, outfile = None, writefits = None )

    # cross ch1 and X-ray
     (power_1x, sig_1x) = cross_power ( amp_ir1, amp, power_ir1, power, f_clip, k_w, writefits = None )
    # cross ch1 ab and X-ray
     (power_1xab, sig_1x) = cross_power ( amp_ab_ir1, amp, power_ab_ir1, power, f_clip, k_w, writefits = None )

    # cross ch2 and X-ray
     (power_2x, sig_2x) = cross_power ( amp_ir2, amp, power_ir2, power, f_clip, k_w, writefits = None )
    # cross ch2 ab and X-ray
     (power_2xab, sig_2x) = cross_power ( amp_ab_ir2, amp, power_ab_ir2, power, f_clip, k_w, writefits = None )

    # cross ch4 and X-ray
     (power_4x, sig_4x) = cross_power ( amp_ir4, amp, power_ir4, power, f_clip, k_w, writefits = None )
    # cross ch4 ab and X-ray
     (power_4xab, sig_4x) = cross_power ( amp_ab_ir4, amp, power_ab_ir4, power, f_clip, k_w, writefits = None )

    # Coherence

def plot_auto_power () :



if __name__ == "__main__":
    main()
