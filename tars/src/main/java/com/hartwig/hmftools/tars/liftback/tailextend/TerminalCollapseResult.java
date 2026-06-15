package com.hartwig.hmftools.tars.liftback.tailextend;

// Result of TerminalMicroJunctionCollapser.tryCollapse. NewStart shifts only on a leading collapse.
public class TerminalCollapseResult
{
    public final boolean Collapsed;
    public final int NewStart;
    public final String NewCigar;

    public TerminalCollapseResult(final boolean collapsed, final int newStart, final String newCigar)
    {
        Collapsed = collapsed;
        NewStart = newStart;
        NewCigar = newCigar;
    }

    public static TerminalCollapseResult unchanged()
    {
        return new TerminalCollapseResult(false, 0, null);
    }
}
