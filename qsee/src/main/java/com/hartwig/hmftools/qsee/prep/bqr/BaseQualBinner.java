package com.hartwig.hmftools.qsee.prep.bqr;

import java.util.EnumMap;

import com.hartwig.hmftools.common.sequencing.SequencingType;

public class BaseQualBinner
{
    private final SequencingType mSequencingType;
    private final EnumMap<BaseQualBin, BaseQualBinRange> mBinRanges;

    public BaseQualBinner(SequencingType sequencingType)
    {
        mSequencingType = sequencingType;
        mBinRanges = getBinRanges();
    }

    public String binNameFrom(byte qual)
    {
        for(BaseQualBin bin : mBinRanges.keySet())
        {
            BaseQualBinRange range = mBinRanges.get(bin);
            if(qual >= range.lowerBound() && qual <= range.upperBound())
            {
                return range.toString();
            }
        }
        throw new RuntimeException("No bin found for base quality: " + qual);
    }

    public SequencingType sequencingType() { return mSequencingType; }
    public EnumMap<BaseQualBin, BaseQualBinRange> binRanges() { return mBinRanges; }

    private EnumMap<BaseQualBin, BaseQualBinRange> getBinRanges()
    {
        EnumMap<BaseQualBin, Byte> lowerBounds = new EnumMap<>(BaseQualBin.class);
        EnumMap<BaseQualBin, Byte> upperBounds = new EnumMap<>(BaseQualBin.class);

        for(byte qual = 0; qual < Byte.MAX_VALUE; qual++)
        {
            if(BaseQualBin.LOW.withinBin(qual, mSequencingType))
                lowerBounds.putIfAbsent(BaseQualBin.LOW, qual);

            if(BaseQualBin.MEDIUM.withinBin(qual, mSequencingType))
                lowerBounds.putIfAbsent(BaseQualBin.MEDIUM, qual);

            if(BaseQualBin.HIGH.withinBin(qual, mSequencingType))
                lowerBounds.putIfAbsent(BaseQualBin.HIGH, qual);
        }

        BaseQualBin previousBin = null;
        for(BaseQualBin currentBin : lowerBounds.keySet())
        {
            if(previousBin != null)
            {
                byte currentLowerBound = lowerBounds.get(currentBin);
                byte previousUpperBound = (byte) (currentLowerBound - 1);
                upperBounds.put(previousBin, previousUpperBound);
            }

            previousBin = currentBin;
        }
        upperBounds.put(BaseQualBin.HIGH, Byte.MAX_VALUE);

        EnumMap<BaseQualBin, BaseQualBinRange> binRanges = new EnumMap<>(BaseQualBin.class);
        for(BaseQualBin bin : lowerBounds.keySet())
        {
            binRanges.put(bin, new BaseQualBinRange(bin, lowerBounds.get(bin), upperBounds.get(bin)));
        }

        return binRanges;
    }
}
