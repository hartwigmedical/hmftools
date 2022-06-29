package com.hartwig.hmftools.svtools.sv_prep;

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
}
