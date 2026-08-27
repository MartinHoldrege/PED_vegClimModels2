# Phenology data notes

Important info from the metadata: 

SOST: start of season time (SOST): NAs: a cell value of 1000 represents water bodies, and a cell value of -1000 represents an area where SOST could not be detected 

SOSN: 255 represents water bodies, and a cell value of 100 represents the area where a SOSN could not be detected due to insufficient change in time-series NDVI, or insufficient input data.To revert to unscaled NDVI, a scale factor [(SOSN-100) * 0.01]

EOST: Valid data ranges from 1 to 450. value of 1000 represents water bodies and a cell value of -1000 represents the area where the EOST could not be detected

EOSN: The EOSN is based on NDVI values (unitless) with a scaling factor applied, and valid values range from 101 to 200. The same scale factor as mentioned in SOSN is also used for EOSN.\

MAXT:  MAXT is day of the year and the values range from 1 to 365

MAXN: The MAXN is based on NDVI values (unitless), and values range from 101 to 200. The same scale factor as mentioned in SOSN is also used for MAXN


DUR:  Note--we don't need to model this, it's derived!

Amplitude (AMP): .[Derived!]

Time Integrated NDVI (TIN): The TIN is unitless and valid values range from 1 to 200. In the TIN data layer, a cell value of 255 represents water bodies and a cell value of 0 represents the area where TIN 





