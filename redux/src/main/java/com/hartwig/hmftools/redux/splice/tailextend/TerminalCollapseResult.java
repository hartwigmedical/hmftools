package com.hartwig.hmftools.redux.splice.tailextend;

// Result of TerminalMicroJunctionCollapser.tryCollapse. Collapsed=false means no change (NewStart and
// NewCigar are unset). NewStart shifts only when the leading side was collapsed.
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
