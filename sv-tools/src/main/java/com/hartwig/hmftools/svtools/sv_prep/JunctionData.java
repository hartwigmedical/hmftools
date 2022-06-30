package com.hartwig.hmftools.svtools.sv_prep;

import java.util.Comparator;

class JunctionData
{
    public final int Position;
    public final byte Orientation;
    public int ExactSupport;
    public int CandidateSupport;

    public JunctionData(final int position, final byte orientation)
    {
        Position = position;
        Orientation = orientation;
        ExactSupport = 1;
        CandidateSupport = 0;
    }

    public int totalSupport() { return ExactSupport + CandidateSupport; }

    public static class JunctionDataSorter implements Comparator<JunctionData>
    {
        // sorts by support descending
        public int compare(final JunctionData first, final JunctionData second)
        {
            if(first.totalSupport() != second.totalSupport())
                return first.totalSupport() < second.totalSupport() ? 1 : -1;

            return 0;
        }
    }

}
