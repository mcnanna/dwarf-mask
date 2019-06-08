#!/usr/bin/env python
import os,sys
from os.path import join,abspath,split
import inspect
from collections import OrderedDict as odict

import numpy as np
from numpy.lib.recfunctions import stack_arrays
import fitsio

import ugali.utils.projector
from ugali.utils.projector import gal2cel, cel2gal
import ugali.utils.idl
from ugali.utils.healpix import ang2pix
from ugali.utils.shell import get_ugali_dir, get_cat_dir
from ugali.utils.logger import logger

import astropy.table
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1


FILTER_DWARFS=True
        

class SourceCatalog(object):
    DATADIR=get_cat_dir()
 
    def __init__(self, filename=None):
        columns = [('name',object),
                   ('ra',float),
                   ('dec',float),
                   ('glon',float),
                   ('glat',float),
                   ('radius', float)] # degrees
        self.data = np.recarray(0,dtype=columns)
        self._load(filename)
        if np.isnan([self.data['glon'],self.data['glat']]).any():
            raise ValueError("Incompatible values")
 
    def __getitem__(self, key):
        """ 
        Support indexing, slicing and direct access.
        """
        try:
            return self.data[key]
        except ValueError as message:
            if key in self.data['name']:
                return self.data[self.data['name'] == key]
            else:
                raise ValueError(message)
 
    def __add__(self, other):
        ret = SourceCatalog()
        ret.data = np.concatenate([self.data,other.data])
        return ret
        
    def __len__(self):
        """ Return the length of the collection.
        """
        return len(self.data)
 
    def _load(self, filename):
        pass
 
    def match(self,lon,lat,coord='gal',tol=0.1,nnearest=1):
        if coord.lower() == 'cel':
            glon, glat = cel2gal(lon,lat)
        else:
            glon,glat = lon, lat
        return ugali.utils.projector.match(glon,glat,self['glon'],self['glat'],tol,nnearest)
   


class McConnachie12(SourceCatalog):
    """
    Catalog of nearby dwarf spheroidal galaxies.
    http://arxiv.org/abs/1204.1562

    http://www.astro.uvic.ca/~alan/Nearby_Dwarf_Database_files/NearbyGalaxies.dat
    """
    def _load(self, filename):
        raw = Vizier.get_catalogs('J/AJ/144/4')[0]

        self.data.resize(len(raw))
        self.data['name'] = raw['Name']

        ra = np.array([map(float, hms.split()) for hms in raw['RAJ2000']])
        dec = np.array([map(float, dms.split()) for dms in raw['DEJ2000']])
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)

        self.data['radius'] = raw['R1']/60.0 # Half-light radius along major axis
        # Could also include elliticity 'Ell' and position angle 'PA'

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat


class McConnachie15(SourceCatalog):
    """
    Catalog of nearby dwarf spheroidal galaxies. Updated September 2015.
    http://arxiv.org/abs/1204.1562

    http://www.astro.uvic.ca/~alan/Nearby_Dwarf_Database_files/NearbyGalaxies.dat
    """
    def _load(self,filename):
        if filename is None: 
            filename = os.path.join(self.DATADIR,"J_AJ_144_4/NearbyGalaxies.dat")
        self.filename = filename
 
        delimiter = [19,3,3,5,3,3,3,6,6,5,5,7,5,5,5,4,4,6,5,5,5,5,5,5,4,4,7,6,6,5,5,5,5,5,5,5,5,5,5,6,5,5,6,5,5,2]
        raw = np.genfromtxt(filename,delimiter=delimiter,usecols=list(range(7))+[26],dtype=['|S19']+7*[float],skip_header=36)

        raw[['LMC' in name for name in raw['f0']].index(True)]['f7'] = 540.0 # LMC radius = 9 deg
        raw[['SMC' in name for name in raw['f0']].index(True)]['f7'] = 180.0 # LMC radius = 3 deg
        raw[['Bootes III' in name for name in raw['f0']].index(True)]['f7'] = 60.0 # Bootes III radius = 1 deg

        self.data.resize(len(raw))
        self.data['name'] = np.char.lstrip(np.char.strip(raw['f0']),'*')

        ra = raw[['f1','f2','f3']].view(float).reshape(len(raw),-1)
        dec = raw[['f4','f5','f6']].view(float).reshape(len(raw),-1)
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)

        radius = raw['f7']
        # McConnachie sets "bad" radius to either 99.99 or 9.999
        bad_r = np.nan # arcmin
        radius = np.array([r if set(str(r))!=set('9.') else bad_r for r in radius]) 
        self.data['radius'] = radius/60.0 # Half-light radius along major axis
        # Could also include ellipticity and position angle

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat


class Rykoff14(SourceCatalog):
    """
    Catalog of red-sequence galaxy clusters.
    http://arxiv.org/abs/1303.3562

    """

    def _load(self, filename):
        if filename is None: 
            filename = os.path.join(self.DATADIR,"redmapper/dr8_run_redmapper_v5.10_lgt20_catalog.fit")
        self.data['ra'] = ugali.utils.projector.hms2dec(raw['RAJ2000'])
        self.filename = filename

        raw = fitsio.read(filename,lower=True)

        self.data.resize(len(raw))
        self.data['name'] = np.char.mod("RedMaPPer %d",raw['mem_match_id'])
        self.data['ra'] = raw['ra']
        self.data['dec'] = raw['dec']
        glon,glat = cel2gal(raw['ra'],raw['dec'])
        self.data['glon'],self.data['glat'] = glon, glat


class Harris96(SourceCatalog):
    """
    Catalog of Milky Way globular clusters.
    Harris, W.E. 1996, AJ, 112, 1487

    http://physwww.physics.mcmaster.ca/~harris/mwgc.dat

    NOTE: There is some inconsistency between Equatorial and
    Galactic coordinates in the catalog. Equatorial seems more
    reliable.
    """
    def _load(self,filename):
        if filename is None: 
            filename = os.path.join(self.DATADIR,"VII_202/mwgc.dat")
        self.filename = filename

        kwargs = dict(delimiter=[12,12,3,3,6,5,3,6,8,8,6],dtype=2*['S12']+7*[float],skip_header=72,skip_footer=363)
        raw = np.genfromtxt(filename,**kwargs)
        
        kwargs2 = dict(delimiter=[13,6,6,8,8,8,9,6,6,8,7,7,5],dtype=['S12']+12*[float], skip_header=433, skip_footer=2)
        raw2 = np.genfromtxt(filename,**kwargs2)

        self.data.resize(len(raw))
        self.data['name'] = np.char.strip(raw['f0'])

        ra = raw[['f2','f3','f4']].view(float).reshape(len(raw),-1)
        dec = raw[['f5','f6','f7']].view(float).reshape(len(raw),-1)

        radius = raw2['f8']/60.0 # Half-light radius
        # Could also use core radius raw2['f7']

        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat

        self.data['radius'] = radius


class Corwen04(SourceCatalog):
    """
    Modern compilation of the New General Catalogue and IC
    """
    def _load(self,filename):
        kwargs = dict(delimiter=[1,1,4,15,3,3,8,3,3,7],usecols=[1,2]+list(range(4,10)),dtype=['S1']+[int]+6*[float])
        if filename is None: 
            raw = []
            for basename in ['VII_239A/ngcpos.dat','VII_239A/icpos.dat']:
                filename = os.path.join(self.DATADIR,basename)
                raw.append(np.genfromtxt(filename,**kwargs))
            raw = np.concatenate(raw)
        else:
            raw = np.genfromtxt(filename,**kwargs)
        self.filename = filename

        #raw = raw[(raw['f0'] != 'I') | (raw['f1'] != 1613)] # IC 1613

        # Some entries are missing...
        raw['f4'] = np.where(np.isnan(raw['f4']),0,raw['f4'])
        raw['f7'] = np.where(np.isnan(raw['f7']),0,raw['f7'])

        self.data.resize(len(raw))
        names = np.where(raw['f0'] == 'N', 'NGC %04i', 'IC %04i')
        self.data['name'] = np.char.mod(names,raw['f1'])

        ra = raw[['f2','f3','f4']].view(float).reshape(len(raw),-1)
        dec = raw[['f5','f6','f7']].view(float).reshape(len(raw),-1)
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat


class Nilson73(SourceCatalog):
    """
    Modern compilation of the Uppsala General Catalog

    http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/26D
    """
    def _load(self,filename):
        raw = Vizier.get_catalogs('VII/26D')[0]
        if FILTER_DWARFS:
            raw = raw[raw['UGC'] != 10822] # Draco (UGC 10822)
            raw = raw[raw['UGC'] != 9749] # Ursa Minor (UGC 9749)
            raw = raw[raw['UGC'] != 5470] # Leo I (UGC 5470)
            raw = raw[raw['UGC'] != 6253] # Leo II (UGC 6253)
        
        self.data.resize(len(raw))
        self.data['name'] = np.char.mod('UGC %s', raw['UGC'])

        ra = np.array([map(float, hms.split()) for hms in raw['RA1950']])
        dec = np.array([map(float, dms.split()) for dms in raw['DE1950']])
        ra = np.vstack([ra.T,np.zeros(len(raw))]).T
        dec = np.vstack([dec.T,np.zeros(len(raw))]).T
        ra1950 = ugali.utils.projector.hms2dec(ra)
        dec1950 = ugali.utils.projector.dms2dec(dec)
        ra2000,dec2000 = ugali.utils.idl.jprecess(ra1950,dec1950)
        self.data['ra'] = ra2000
        self.data['dec'] = dec2000
        
        self.data['radius'] = raw['MajAxis']/60.0 # Major axis on POSS blue print
        # Could also use minor axis 'MinAxis' (no position angle given)

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat

class Webbink85(SourceCatalog):
    """
    Structure parameters of Galactic globular clusters
    http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/151

    NOTE: Includes Reticulum and some open clusters
    http://spider.seds.org/spider/MWGC/mwgc.html
    """
    def _load(self,filename):
        catalog = Vizier.get_catalogs('VII/151')
        raw = astropy.table.vstack([catalog[0], catalog[2]], metadata_conflicts='silent')
      
        self.data.resize(len(raw))
        self.data['name'] = np.char.join(' ', np.char.split(raw['Name1']))

        ra = np.array([map(float, hms.split()) for hms in raw['RA1950']])
        dec = np.array([map(float, dms.split()) for dms in raw['DE1950']])
        dec = np.vstack([dec.T,np.zeros(len(raw))]).T
        ra1950 = ugali.utils.projector.hms2dec(ra)
        dec1950 = ugali.utils.projector.dms2dec(dec)
        ra2000,dec2000 = ugali.utils.idl.jprecess(ra1950,dec1950)
        self.data['ra'] = ra2000
        self.data['dec'] = dec2000
        
        self.data['radius'] = (10**np.array(raw['lgtt']))/60.0 # log(radius)
        # Could also use log(core radius) 'lgtc'

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat


class Kharchenko13(SourceCatalog):
    """
    Global survey of star clusters in the Milky Way
    http://vizier.cfa.harvard.edu/viz-bin/Cat?J/A%2bA/558/A53

    NOTE: CEL and GAL coordinates are consistent to < 0.01 deg.
    """
    def _load(self,filename):
        raw = Vizier.get_catalogs('J/A+A/558/A53')
        if FILTER_DWARFS:
            raw = [ raw[0][raw[0]['MWSC'] != '2020'] , raw[1], raw[2] ] # Coma Berenices (Melotte_111)

        self.data.resize(len(raw[0]))

        ids = np.array(map(int, raw[0]['MWSC'])) - 1
        self.data['name'] = raw[1]['Name'][ids]

        self.data['ra'] = raw[0]['RAJ2000']
        self.data['dec'] = raw[0]['DEJ2000']

        self.data['radius'] = raw[0]['r1'] # Radius of central part
        # Could also use core radius 'r0' or radius 'r2'

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat



class Bica08(SourceCatalog):
    """
    LMC star clusters
    http://cdsarc.u-strasbg.fr/viz-bin/Cat?J/MNRAS/389/678

    NOTE: CEL and GAL coordinates are consistent to < 0.01 deg.
    """
    def _load(self,filename):
        raw = Vizier.get_catalogs('J/MNRAS/389/678')[0]

        self.data.resize(len(raw))
        self.data['name'] = np.char.strip(raw['Names'])

        ra = np.array([map(float, hms.split()) for hms in raw['RAJ2000']])
        dec = np.array([map(float, dms.split()) for dms in raw['DEJ2000']])
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)
        
        self.data['radius'] = raw['amaj']/60.0 # Major axis
        # Could also use minor axis 'amin' and position angle 'PA'

        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat

class WEBDA14(SourceCatalog):
    """
    Open cluster database.
    http://www.univie.ac.at/webda/cgi-bin/selname.cgi?auth=
    
    """
    def _load(self,filename):
        kwargs = dict(delimiter='\t',usecols=[0,1,2,9],dtype=['S18',float,float,float])
        if filename is None: 
            filename = os.path.join(self.DATADIR,"WEBDA/webda.tsv")
        self.filename = filename
        raw = np.genfromtxt(filename,**kwargs)
        
        self.data.resize(len(raw))
        self.data['name'] = np.char.strip(raw['f0'])

        self.data['glon'] = raw['f1']
        self.data['glat'] = raw['f2']

        self.data['radius'] = raw['f3']/2.0/60.0 # Diameter

        ra,dec = gal2cel(self.data['glon'],self.data['glat'])
        self.data['ra'],self.data['dec'] = ra,dec

class ExtraDwarfs(SourceCatalog):
    """
    Collection of dwarf galaxy candidates discovered in 2015
    """
    def _load(self,filename):
        kwargs = dict(delimiter=',')
        if filename is None: 
            filename = os.path.join(self.DATADIR,"extras/extra_dwarfs.csv")
        self.filename = filename
        raw = np.recfromcsv(filename,**kwargs)
        
        self.data.resize(len(raw))
        self.data['name'] = raw['name']
        
        self.data['ra'] = raw['ra']
        self.data['dec'] = raw['dec']

        self.data['radius'] = raw['rhalf']/60.0 # Half-light radius, arcmin

        self.data['glon'],self.data['glat'] = cel2gal(raw['ra'],raw['dec'])

class ExtraClusters(SourceCatalog):
    """
    Collection of recently discovered star clusters
    """
    def _load(self,filename):
        kwargs = dict(delimiter=',')
        if filename is None: 
            filename = os.path.join(self.DATADIR,"extras/extra_clusters.csv")
        self.filename = filename
        raw = np.recfromcsv(filename,**kwargs)
        if FILTER_DWARFS:
            raw = raw[raw['name'] != 'Kim 2'] # Indus 1 (Kim 2)
            pass
        
        self.data.resize(len(raw))
        self.data['name'] = raw['name']
        
        self.data['ra'] = raw['ra']
        self.data['dec'] = raw['dec']

        self.data['radius'] = raw['rhalf']/60.0 # Half-light radius, arcmin

        self.data['glon'],self.data['glat'] = cel2gal(raw['ra'],raw['dec'])

class ExtraStructures(SourceCatalog):
    """
    Collection of recently discovered structures
    """
    def _load(self,filename):
        kwargs = dict(delimiter=',')
        if filename is None: 
            filename = os.path.join(self.DATADIR,"extras/extra_structures.csv")
        self.filename = filename
        raw = np.recfromcsv(filename,**kwargs)
        
        self.data.resize(len(raw))
        self.data['name'] = raw['name']
        
        self.data['ra'] = raw['ra']
        self.data['dec'] = raw['dec']

        self.data['radius'] = raw['width'] # Only applies to ATLAS stream

        self.data['glon'],self.data['glat'] = cel2gal(raw['ra'],raw['dec'])


def catalogFactory(name, **kwargs):
    """
    Factory for various catalogs.
    """
    fn = lambda member: inspect.isclass(member) and member.__module__==__name__
    catalogs = odict(inspect.getmembers(sys.modules[__name__], fn))

    if name not in list(catalogs.keys()):
        msg = "%s not found in catalogs:\n %s"%(name,list(kernels.keys()))
        logger.error(msg)
        msg = "Unrecognized catalog: %s"%name
        raise Exception(msg)

    return catalogs[name](**kwargs)

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('args',nargs=argparse.REMAINDER)
    opts = parser.parse_args(); args = opts.args

