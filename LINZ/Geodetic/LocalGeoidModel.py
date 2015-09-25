# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from .Adjustment import Plugin, Options

def LocalGeoidModelPlugin( Plugin ):
        
    def getOptions():
        return Options(localGeoidModel=None)

    def setConfigOption( self, item, value ):
        if item == 'local_geoid':
            try:
                gparts=value.split()
                code=gparts.pop(0)
                geoidHeight,xi,eta=(float(v) for v in gparts[0:3])
                range=float(gparts[4]) if len(gparts)==5 else None
            except:
                raise RuntimeError("Invalid local_geoid: "+value)
            self.options.localGeoidModel={'code':code,'height':geoidHeight,'xi':xi,'eta':eta, 'range':range}
            return True
        return False

    def setupLocalGeoid( self, geoidModel ):
        refStationCode=geoidModel['code']
        refStation=self.stations.get( geoidModel['code'] )
        if refStation is None:
            raise RuntimeError("Local geoid model reference station {0} is not defined"
                               .format(refStationCode))
        gheight=geoidModel['height']
        xi=geoidModel['xi']/3600.0
        eta=geoidModel['eta']/3600.0
        xirad=math.radians(xi)
        etarad=math.radians(eta)
        range=geoidModel.get('range')
        if range is not None:
            range = range*range

        xyz0=refStation.xyz()
        enu_axes=refStation.enu()

        for s in self.network.stations:
            denu=enu.dot(s.xyz()-xyz0)
            if range is not None and (denu[0]*denu[0]+denu[1]*denu[1]) > range:
                continue
            s.setXiEta([xi,eta])

    def preSetupParameters()
        # Apply a geoid model
        if options.localGeoidModel is not None:
            gm=options.localGeoidModel
            code=gm['code']
            geoidHeight=gm['height']
            xieta=[gm['xi']/3600.0,gm['eta']/3600.0]
            range=gm['range']
            self.stations.setLocalGeoid(code,geoidHeight,xieta,range)
