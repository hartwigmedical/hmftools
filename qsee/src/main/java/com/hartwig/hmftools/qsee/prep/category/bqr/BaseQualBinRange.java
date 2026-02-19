package com.hartwig.hmftools.qsee.prep.category.bqr;

public class BaseQualBinRange
{
    private final BaseQualBin Bin;
    private final byte LowerBound;
    private final byte UpperBound;

    public BaseQualBinRange(BaseQualBin bin, byte lowerBound, byte upperBound)
    {
        Bin = bin;
        LowerBound = lowerBound;
        UpperBound = upperBound;
    }

    public BaseQualBin bin(){ return Bin; }
    public byte lowerBound(){ return LowerBound; }
    public byte upperBound(){ return UpperBound; }

    @Override
    public String toString()
    {
        return (Bin == BaseQualBin.HIGH)
                ? String.format("%s (%d+)", Bin, LowerBound)
                : String.format("%s (%d-%d)", Bin, LowerBound, UpperBound);
    }
}
