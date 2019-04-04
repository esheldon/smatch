import os
import numpy as np
import biggles
import esutil as eu
import argparse
import smatch
import time

parser=argparse.ArgumentParser()
parser.add_argument('--ntrial',type=int,default=100)
parser.add_argument('--density-arcmin',type=float,default=30)
parser.add_argument('--area-arcmin',type=float,default=40)

parser.add_argument('--seed',type=int,default=None)
parser.add_argument('--maxmatch',type=int,default=1)

parser.add_argument('--show',action='store_true')

def get_data(args):

    area_deg = args.area_arcmin/3600.0

    linear_size = np.sqrt(area_deg)

    num = int( args.area_arcmin*args.density_arcmin )

    return eu.coords.randsphere(
        num,
        ra_range=[100,100+linear_size],
        dec_range=[-0.5*linear_size,0.5*linear_size],
    )


def match_smatch(args,
                 nside,
                 ra1, dec1, radius1, ra2, dec2):

    cat = smatch.Catalog(ra1, dec1, radius1, nside=nside)
    cat.match(ra2, dec2, maxmatch=args.maxmatch)

def doplot(plt, args, nsides, times, radius_arcsec, color):

    #times = np.array(times)
    #times /= times.max()

    label='rad %g arcsec' % radius_arcsec
    plt.add(
        biggles.Curve(nsides, times, label=label, color=color),
    )

def main():
    args=parser.parse_args()
    print('seed:',args.seed)

    np.random.seed(args.seed)


    #nsides=[512,1024,2048,4096,8192]
    nsides=np.arange(500,8200,500)

    key=biggles.PlotKey(0.1, 0.9, halign='left')
    plt=biggles.FramedPlot(
        xlabel='nside',
        ylabel='time [s]',
        aspect_ratio=1,
        xrange=[0.0, 1.05*nsides[-1]],
        #yrange=[0.002,0.1],
        yrange=[0.0,0.04],
        #ylog=True,
        key=key,
    )
    radii_arcsec = [
        2.0,
        10.0,
        60.0,
        200.0,
    ]
    colors=[
        'black',
        'darkgreen',
        'blue',
        'red',
    ]


    for radius_arcsec,color in zip(radii_arcsec,colors):
        print('radius:',radius_arcsec)
        radius1 = radius_arcsec/3600.0

        times=[]
        for nside in nsides:
            ttimes = []
            for i in range(args.ntrial):
                ra1, dec1 = get_data(args)
                ra2, dec2 = get_data(args)

                tm0 = time.time()
                match_smatch(args, nside, ra1, dec1, radius1, ra2, dec2)
                tm = time.time()-tm0
                ttimes.append(tm)

            tm = min(ttimes)

            times.append(tm)

        doplot(plt, args, nsides, times, radius_arcsec, color)

    fname='smatch-times-density-%g.png' % args.density_arcmin
    print('writing:',fname)
    plt.write(fname, dpi=150)

main()
