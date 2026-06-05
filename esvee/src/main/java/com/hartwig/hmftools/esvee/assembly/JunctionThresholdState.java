package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

public class JunctionThresholdState
{
    public boolean ExtensionLengthValid;
    public boolean SecondExtensionLengthValid;
    public boolean MinReadsValid;
    public boolean UsesLowerSagaLimits;
    public boolean IsLINE;
    public boolean RequireNonExtensionSupport;

    // thresholds
    public int MinRequiredReads;
    public int MinExtensionLength;
    public int MinSecondExtensionLength;

    public JunctionThresholdState()
    {
        ExtensionLengthValid = false;
        SecondExtensionLengthValid = false;
        MinReadsValid = false;
        IsLINE = false;
        RequireNonExtensionSupport = false;
        MinRequiredReads = 0;
        MinExtensionLength = 0;
        MinSecondExtensionLength = 0;
        UsesLowerSagaLimits = false;
    }

    public boolean isValid() { return ExtensionLengthValid && SecondExtensionLengthValid && MinReadsValid; }

    public String toString()
    {
        return format("minReads(%s req=%d) mainLength(%s req=%d) secondLength(%s req=%d) line(%s) allowJS(%s) sagaLimits(%s)",
                MinReadsValid, MinRequiredReads, ExtensionLengthValid, MinExtensionLength,
                SecondExtensionLengthValid, MinSecondExtensionLength, IsLINE, RequireNonExtensionSupport, UsesLowerSagaLimits);
    }
}
