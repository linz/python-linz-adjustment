# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from .Adjustment import Plugin
from . import StationLocator

class StationLocatorPlugin( Plugin ):

    def pluginOptions( self ):
        return dict(
            calcMissingCoords=False,
            debugCalcMissingCoords=False
        );

    def setConfigOption( self, item, value ):
        if item == 'calculate_missing_stations':
            self.options.calcMissingCoords=self.options.boolOption(value)
        elif item == 'debug_calculate_missing_stations':
            self.options.debugCalcMissingCoords=self.options.boolOption(value)
        else:
            return False
        return True

    def preSetup( self ):
        options=self.options
        # Calculate missing stations
        if options.calcMissingCoords:
            from . import StationLocator
            adjustment=self.adjustment
            write=None
            if options.debugCalcMissingCoords:
                write=adjustment.write
                write("\nCalculating missing station coordinates\n")
            nupdated=StationLocator.locateStations(
                adjustment.stations,
                adjustment.observations,
                write)
            if nupdated > 0:
                adjustment.write("\nApproximate coordinates calculated for {0} stations\n"
                           .format(nupdated))
