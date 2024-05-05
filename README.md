# LINZ.adjustment package

The python-linz-adjustment software is used to calculate least squares adjustment
of survey observations for quality control of the observations and to calculate
coordinates of the observed survey stations.  This module uses the linear algebra
functions of the python numpy module to set up and solve the least squares equations.

While this is a general purpose least squares adjustment module it was written primarily
to support the `python-linz-pyaxis <https://github.com/linz/python-linz-pyaxis>`_ module
used to analyse data for ITRF local tie surveys.

The software uses a plugin architecture that allows additional funtionality to be added.
For example a plugin can add additional parameters to the adjustment.  The pyaxis software
is implemented as a plugin to this software.

The adjustment software is designed for relatively small adjustments - it readily handles
adjustments of a few hundred parameters.

## Running the software

The software uses three types of input file:

* observation file.  The adjustment can include multiple observation files.  Each file defines survey observations.  The software supports a number of data formats including a CSV (comma-separated value) format,
`SINEX <https://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html>`_ data of station coordinates from GNSS processing.  Additional observation formats can be supported using plugins.
* coordinate file.  Defines the stations in the adjustment and their approximate coordinates.  This is not required if the software can calculate all the station coordinates from the input data using the station_locator plugin.  
of any stations defined in the input data which for which an approximate coordinate is not defined
* configuration file.  Defines the configuration of the adjustment, including the names of the observation and coordinate files and the names of output files created by the adjustment.  The software can create a sample configuration which can be edited and used to run the adjustment.

The software is run using the snappy command.  To create an example configuration file use:

  snappy -c myconfig.cfg
  
This will create a sample configuration file called myconfig.cfg.  This can be edited and then used to run an adjustment using the command:

  snappy myconfig.cfg
  
## Observation files

In order to run an adjustment the first thing to do is create input files of the survey observation data.  Two data file formats are supported:

* CSV format.  Can define horizontal angle, horizontal distance, zenith distance, slope distance, azimuth, and levelling observations.  Details are below
* SINEX format.  Defines GNSS XYZ coordinate observations.

Observation file are defined in the configuration file by data_file configuration commmands, for example

```text
   data_file control_survey.csv
```

The file type is defined by the extension (csv or snx).  If the file has a different extension the filename can be preceded by the type eg `csv:control_survey.dat` to define a CSV formatted file named control_survey.dat.  The data_file command can also define additional parameters of the data format using param=value options, where param is the name of the parameter and value is its value.  The available parameters for each format are defined below.

Each station in the adjustment is identified by a code, eg `WARK1` which is unique in the survey.  This is used to reference the station in the input and output files.

### CSV file format

The CSV format is used to define distance, angle, and levelling survey observations. The observation types supported are:

* `HA`: horizontal angle observations in decimal degrees
* `AZ`: azimuth observations in decimal degrees (not usually observed - may be included as dummy observations to constrain an adjustment)
* `SD`: slope distance observations in metres
* `HD`: horizontal distance observations in metres (not usually observed - may be output by total stations)
* `ZD`: zenith distance observations measured in decimal degrees from the vertical axis
* `LV`: levelling height difference observations in metres

The first line of the file defines the names of each field in the file.  The expected field names are:

* `fromstn` - the code of the survey station from which the file is made
* `fromhgt` - the height of the equipment above the station
* `tostn` - the code of the survey station to which the observation is made
* `tohgt` - the height of the target above the station
* `date` - the date and time of the observation in format `yyyy-mm-dd hh:mm:ss`
* `obsset` - for observations made in a set (horizontal angles) this identifies the set to which the observation belongs
* `value` - the observation value
* `error` - the uncertainty of the observation

As an alternative to specifying the observation type the CSV file can contain columns named by the data type, for example
`ha_value,ha_error,sd_value,sd_error`.

Additional fields can be read from the data file can be read by adding `attributes=fieldname1+fieldname2+...` in the data_file command.  Some adjustment options can use these attribute values (for example in the setup_heights plugin).  They can also be output to the residuals CSV file, which may be useful for analysis of the adjustment results.

## Coordinate file

The configuration file may define one or more coordinate files which define the initial (trial) locations of stations, or the
fixed coordinates of reference stations. Coordinates are defined either as geocentrix X,Y,Z coordinates or longitude,latitude,ellipsoidal height coordinates.  Longitude and latitude values are in decimal degrees.

The coordinate file is a CSV file with the following columns:

* `code` - the code used to identify the stations in the data files and adjustment
* `name` - (optional) a descriptive name for the station
* `x` or `lon` or `longitude` - the X or latitude ordinate of the station
* `y` or `lat` or `latitude` - the Y or latitude ordinate of the station
* `z` or `hgt` or `height` or `ellheight` - the Z or ellipsoidal height ordinate of the station
* `geoid` - (optional) the geoid height at the station in metres
* `xi` - (optional) the north-south deflection of the vertical in seconds
* `eta` - (optional) the east-west deflection of the vertical in seconds

## Plugins

The following plugins can be used in the adjustment.  The plugin can be included in the configuration file using a command such as:

  use_plugin local_geoid_model
  
Plugins can also be included by adding `-p plugin_name1+plugin_name2..` to the snappy command line.  To create a sample configuration file including the plugin configuration parameters run snappy with the plugin options on the command line and the -c option to create an example configuration file.

### Local geoid model plugin

Defines a local geoid model.  Defines a simplistic geoid model defined by a geoid height and deflection of the vertical at a reference station.  The geoid height is interpolated from the reference station based on the deflection.  This is only appropriate for a survey of limited extents.

### Setup height plugin

Adds functionality to calculate the instrument/target height for each setup.  Requires attributes are loaded from the CSV file identifying the instrument/target setup for each station.  These are then used to identify the stations for which a setup is calculated.  By default the identified setup heights will be calculated if this plugin is enabled.  Specific setups can be fixed to a specified value, or set to a measured value with an uncertainty

### Station locator plugin

Adds functionality to calculate the coordinates of stations not defined in the coordinate file.

## Build notes

NOTE: Currently debian package fails on testing as doesn't find LINZ.Geodetic.Ellipsoid.  Build with `DEB_BUILD_OPTIONS=nocheck`.
