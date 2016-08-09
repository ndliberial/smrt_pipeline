merSize               = 14
merThreshold          = 0
merDistinct           = 0.9995
merTotal              = 0.995

doOBT                 = 0
doExtendClearRanges   = 0

unitigger             = bogart

ovlErrorRate          = 0.35  #  Compute overlaps up to 35% error
utgGraphErrorRate     = 0.35  #  Unitigs at 35% error
utgMergeErrorRate     = 0.35  #  Unitigs at 35% error
cnsErrorRate          = 0.35  #  Needed to allow ovlErrorRate=0.35
cgwErrorRate          = 0.35  #  Needed to allow ovlErrorRate=0.35

ovlConcurrency        = 1
cnsConcurrency        = 8

ovlThreads            = 8
ovlHashBits           = 22
ovlHashBlockLength    = 10000000
ovlRefBlockSize       = 25000

#cnsReduceUnitigs      = 0 0   #  Always use only uncontained reads for consensus
cnsReuseUnitigs       = 1     #  With no mates, no need to redo consensus

cnsMinFrags           = 1000
cnsPartitions         = 256
