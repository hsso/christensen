from os.path import expanduser, join

basedir = expanduser('~/HssO/Christensen/')
datadir = join(basedir, 'data')
paperdir = join(basedir, 'paper')
figsdir = join(basedir, 'paper', 'figures')

freq0 = {'H2O': 556.9359877, 'NH3': 572.4980678}

hifi_obsid = 1342204014
obsids = (1342186621, 1342186622, 1342203478, 1342203479, 1342204014)

horizons_file = {1342186621: "horizons_2009.txt", 1342203478: "horizons_2010.txt"}

# ra and dec for 1342186621
radec1 = (287.58229, -9.10377)
radec2 = (287.58256, -9.10535)

radec = (287.58242, -9.10456)

psfdir = expanduser('~/project/CH_Cyg/PSF')
psf_vesta = {'blue': 'psf20_blu_10_vesta_od160_ama+63_recentered.fits',
            'red': 'psf20_red_10_vesta_od160_ama+63_recentered.fits'}
