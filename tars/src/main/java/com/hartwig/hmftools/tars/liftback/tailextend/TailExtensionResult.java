package com.hartwig.hmftools.tars.liftback.tailextend;

// Result of TerminalReconciler.tryExtend. newStart shifts only on a leading extension.
public record TailExtensionResult(
        boolean extended, int newStart, String newCigar, int basesExtendedLead, int basesExtendedTrail)
{
    public static TailExtensionResult unchanged()
    {
        return new TailExtensionResult(false, 0, null, 0, 0);
    }
}
