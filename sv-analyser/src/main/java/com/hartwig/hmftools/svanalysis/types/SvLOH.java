package com.hartwig.hmftools.svanalysis.types;

public class SvLOH
{
    public final String SampleId;
    public final String Chromosome;
    public final int CnIdStart;
    public final int CnIdEnd;
    public final long PosStart;
    public final long PosEnd;
    public final String SegStart;
    public final String SegEnd;
    public final double PrevCN;
    public final double StartCN;
    public final double EndCN;
    public final double MinCN;
    public final int SegCount;
    public final long Length;
    public final String StartSV;
    public final String EndSV;

    public SvLOH(
            final String sampleId, final String chr, int cnIdStart,
            final int cnIdEnd,
            final long posStart,
            final long posEnd,
            final String segStart,
            final String segEnd,
            final double prevCN,
            final double startCN,
            final double endCN,
            final double minCN,
            final int segCount,
            final long length,
            final String startSV,
            final String endSV)
    {
        SampleId = sampleId;
        Chromosome = chr;
        CnIdStart = cnIdStart;
        CnIdEnd = cnIdEnd;
        PosStart = posStart;
        PosEnd = posEnd;
        SegStart = segStart;
        SegEnd = segEnd;
        PrevCN = prevCN;
        StartCN = startCN;
        EndCN = endCN;
        MinCN = minCN;
        SegCount = segCount;
        Length = length;
        StartSV = startSV;
        EndSV = endSV;
    }

}
