# Handmade by DCP
FIELDNAME CygnusA
N_SCANS   1 
N_INPUTS  512
N_CHANS   109           # number of channels in spectrum
CORRTYPE  B             # correlation type to use. 'C'(cross), 'B'(both), or 'A'(auto)
INT_TIME  8.533         # integration time of scan in seconds
FREQCENT  57.468        # observing center freq in MHz
BANDWIDTH 2.616         # total bandwidth in MHz
# To phase to the zenith, these mush be the HA, RA and Dec of the zenith.
HA_HRS    2.13          # NORMALLY 0 the RA of the desired phase centre (hours)
RA_HRS    19.98         # the RA of the desired phase centre (hours) CygA
DEC_DEGS  40.73         # the DEC of the desired phase centre (degs) CygA
DATE      20130822      # YYYYMMDD
TIME      034127        # HHMMSS
INVERT_FREQ 0           # 1 if the freq decreases with channel number
CONJUGATE   0           # conjugate the raw data to fix sign convention problem if necessary
