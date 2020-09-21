from scipy import complex_
import mmap
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import array
from struct import unpack
import PPR_Keywords as PPR

def get_dataformat(dt):
    if dt is 0:
        return 'B'
        # datasize = 1
    elif dt is 1:
        return 'b'
        # datasize = 1
    elif dt is 2:
        return 'h'
        # datasize = 2
    elif dt is 3:
        return 'h'
        # datasize = 2
    elif dt is 4:
        return 'l'
        # datasize = 4
    elif dt is 5:
        return 'f'
        # datasize = 4
    elif dt is 6:
        return 'd'
        # datasize = 8
    else:
        return 'i'
        # datasize = 4
    
def get_mrd_3d(filename_):
    # Open file as read-binary only
    fidf = open(filename_, 'r+b')  # open(filename_,'rb')
    fid = mmap.mmap(fidf.fileno(), 0)
    
    # Read first 4 values from header
    val = unpack('llll', fid.read(16))
    
    # Get dimensions from this
    no_samples, no_views, no_views_2, no_slices = val[0], val[1], val[2], val[3]
    
    # Work out datatype of data
    fid.seek(18)
    datatype_ = unpack('h', fid.read(2))[0]
    fid.seek(48)
    scaling = unpack('f', fid.read(4))
    bitsperpixel = unpack('f', fid.read(4))
    
    fid.seek(152)
    val = unpack('ii', fid.read(8))
    no_echoes = val[0]
    no_expts = val[1]

    fid.seek(256)
    text = fid.read(256)

    dim = [no_expts, no_echoes, no_slices, no_views_2, no_views, no_samples]

    dt = datatype_
   
    if dt >= 0x10:
        iscomplex = 2
        dt = dt-0x10
    else:
        iscomplex = 1
    
    dataformat = get_dataformat(dt)

    # Compute the number of values expected to be read from the dimensions
    num2read = no_expts*no_echoes*no_slices*no_views*no_views_2 * no_samples*iscomplex
   
    fid.seek(512)
    
    m_total = array.array(dataformat)
    m_total.fromfile(fid, num2read)
    if len(m_total) != num2read:
        print("We have a problem...(file length/read mismatch)")
        return 0
    
    par = PPR.ParseKeywords(fid)
    fid.close()

    if iscomplex is 2:
        m_real = m_total[::2]
        m_imag = m_total[1::2]
        m_C = np.vectorize(complex)(m_real, m_imag)
        m_real = None
        m_imag = None
        m_total = None
    else:
        m_C = m_total
        m_total = None

    n = 0
    
    ord_ = list(range(no_views))

    if 'VAR_centric_on' in par:
        if int(par['VAR centric_on']) == 1:
            val = int(no_views * 0.5)
            ord_ = list(range(no_views))
            for x in range(val):
                ord_[2*x] = val + x
                ord_[2*x+1] = val - x - 1
    elif 'VAR pe1_order' in par:
        if int(par['VAR pe1_order']) == 1:
            val = int(no_views * 0.5)
            ord_ = list(range(no_views))
            for x in range(val):
                ord_[2*x] = val + x
                ord_[2*x+1] = val - x - 1
    
    ord2 = list(range(no_views_2))
    
    if 'VAR pe2_centric_on' in par:

        if int(par['VAR pe2_centric_on']) == 1:
            val = int(no_views_2 * 0.5)
            ord2 = list(range(no_views_2))
            for x in range(val):
                ord2[2*x] = val + x
                ord2[2*x + 1] = val - x - 1
            
    im = np.zeros(shape=(no_expts, no_echoes, no_slices, no_views, no_views_2,
                         no_samples), dtype=complex_)

    for a in range(no_expts):
        for b in range(no_echoes):
            for c in range(no_slices):
                for d in range(no_views):
                    for e in range(no_views_2):
                        for f in range(no_samples):
                            im[a][b][c][ord_[d]][ord2[e]][f] = m_C[f+n]
                        n += no_samples

    ord_ = None
    
    return im, dim, par

data = get_mrd_3d(filename_="D:/Dropbox\MitsakosPhD/2D_Rad_FINAL_1_TE5.MRD")
mat = np.squeeze(data[0])
df = pd.DataFrame(mat)
d0 = mat[0,:]
print(d0.shape)
plt.plot(d0)
plt.show()