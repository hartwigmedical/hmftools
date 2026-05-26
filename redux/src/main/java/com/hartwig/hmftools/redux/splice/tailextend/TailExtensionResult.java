package com.hartwig.hmftools.redux.splice.tailextend;

// Result of SoftclipTailExtender.tryExtend. Extended=false means no change (NewStart and NewCigar
// are unset). NewStart shifts only when the leading side was extended.
public class TailExtensionResult
{
    public final boolean Extended;
    public final int NewStart;
    public final String NewCigar;
    public final int BasesExtendedLead;
    public final int BasesExtendedTrail;

    public TailExtensionResult(
            final boolean extended, final int newStart, final String newCigar,
            final int basesExtendedLead, final int basesExtendedTrail)
    {
        Extended = extended;
        NewStart = newStart;
        NewCigar = newCigar;
        BasesExtendedLead = basesExtendedLead;
        BasesExtendedTrail = basesExtendedTrail;
    }

    public static TailExtensionResult unchanged()
    {
        return new TailExtensionResult(false, 0, null, 0, 0);
    }
}
