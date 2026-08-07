package com.hartwig.hmftools.purple.reseg;

public final class PeakTroughData
{
    // a validated peak or trough
    public final double Level;
    public final double Value;
    public final double SupportBelowLevel;
    public final double SupportAboveLevel;

    public PeakTroughData(final double level, final double value, final double supportBelowLevel, final double supportAboveLevel)
    {
        Level = level;
        Value = value;
        SupportBelowLevel = supportBelowLevel;
        SupportAboveLevel = supportAboveLevel;
    }

    public String toString()
    {
        return String.format("level(%.2f) value(%.1f) support(below=%.2f above=%.2f)",
                Level, Value, SupportBelowLevel, SupportAboveLevel);
    }
}
