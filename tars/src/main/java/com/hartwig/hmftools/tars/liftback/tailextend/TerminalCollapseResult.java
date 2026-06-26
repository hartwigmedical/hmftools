package com.hartwig.hmftools.tars.liftback.tailextend;

// Result of TerminalReconciler.tryCollapse. newStart shifts only on a leading collapse.
public record TerminalCollapseResult(boolean collapsed, int newStart, String newCigar)
{
    public static TerminalCollapseResult unchanged()
    {
        return new TerminalCollapseResult(false, 0, null);
    }
}
