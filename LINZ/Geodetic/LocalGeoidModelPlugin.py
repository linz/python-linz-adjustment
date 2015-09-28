# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from .Adjustment import Plugin, Options

class LocalGeoidModelPlugin( Plugin ):
        
    def pluginOptions(self):
        return dict(localGeoidModel=None)

    def setConfigOption( self, item, value ):
        if item == 'local_geoid':
            try:
                gparts=value.split()
                code=gparts.pop(0)
                geoidHeight,xi,eta=(float(v) for v in gparts[0:3])
                range=float(gparts[3]) if len(gparts) > 3 else None
            except:
                raise RuntimeError("Invalid local_geoid: "+value)
            self.options.localGeoidModel={'code':code,'height':geoidHeight,'xi':xi,'eta':eta, 'range':range}
            return True
        return False

    def setupLocalGeoid( self, geoidModel ):
        stations=self.adjustment.stations
        refStationCode=geoidModel['code']
        refStation=stations.get( geoidModel['code'] )
        if refStation is None:
            raise RuntimeError("Local geoid model reference station {0} is not defined"
                               .format(refStationCode))
        gheight=geoidModel['height']
        xieta=np.array([geoidModel['xi']/3600.0, geoidModel['eta']/3600.0])
        slope=-np.radians(xieta)

        range=geoidModel.get('range')
        if range is not None:
            range = range*range

        xyz0=refStation.xyz()
        enu=refStation.enu()

        # Note slope is dh/dn, dh/de as xieta are ordered lat,lon

        for s in stations.stations():
            denu=enu.dot(s.xyz()-xyz0)
            if range is not None and (denu[0]*denu[0]+denu[1]*denu[1]) > range:
                continue
            ghgt=gheight+slope.dot((denu[1],denu[0]))
            s.setXiEta(xieta)
            s.setGeoidHeight(ghgt)

    def preSetupParameters( self ):
        # Apply a geoid model
        if self.options.localGeoidModel is not None:
            self.setupLocalGeoid(self.options.localGeoidModel)
