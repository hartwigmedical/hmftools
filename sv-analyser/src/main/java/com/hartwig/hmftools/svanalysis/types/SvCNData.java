package com.hartwig.hmftools.svanalysis.types;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;

public class SvCNData {

    final private int mId;
    final public String Chromosome;
    final public long StartPos;
    final public long EndPos;
    final public double CopyNumber;
    final public String SegStart;
    final public String SegEnd;
    final public String Method;
    final public int BafCount;
    final public double ObservedBaf;
    final public double ActualBaf;
    final public int DepthWindowCount;

    private int mIndex; // in the source table

    private StructuralVariantData mSvData; // linked if known
    private boolean mSvLinkOnStart;

    public SvCNData(final int Id, final String Chromosome, long StartPos, long EndPos,
            final String SegStart, final String SegEnd, int bafCount, double observedBaf, double actualBaf,
            double copyNumber, final String method)
    {
        mId = Id;
        this.Chromosome = Chromosome;
        this.StartPos = StartPos;
        this.EndPos = EndPos;
        CopyNumber = copyNumber;
        this.SegStart = SegStart;
        this.SegEnd = SegEnd;
        Method = method;
        BafCount = bafCount;
        ObservedBaf = observedBaf;
        ActualBaf = actualBaf;
        mSvData = null;
        mSvLinkOnStart = false;
        DepthWindowCount = 0;
        mIndex = 0;
    }

    public SvCNData(final PurpleCopyNumber record, int id)
    {
        mId = id;
        Chromosome = record.chromosome();
        StartPos = record.start();
        EndPos = record.end();
        CopyNumber = record.averageTumorCopyNumber();
        SegStart = record.segmentStartSupport().toString();
        SegEnd = record.segmentEndSupport().toString();
        Method = record.method().toString();
        BafCount = record.bafCount();
        ObservedBaf = record.averageObservedBAF();
        ActualBaf = record.averageActualBAF();
        DepthWindowCount = record.depthWindowCount();
        mIndex = 0;
    }

    public int id() { return mId; }
    public final String chromosome() { return Chromosome; }
    public long startPos() { return StartPos; }
    public long endPos() { return EndPos; }
    public long position(boolean useStart) { return useStart ? StartPos : EndPos; }
    public double copyNumber() { return CopyNumber; }
    public final String segStart() { return SegStart; }
    public final String segEnd() { return SegEnd; }
    public final String method() { return Method; }
    public int bafCount() { return BafCount; }
    public double observedBaf() { return ObservedBaf; }
    public double actualBaf() { return ActualBaf; }

    public int getIndex() { return mIndex; }
    public void setIndex(int index) { mIndex = index; }

    public void setStructuralVariantData(final StructuralVariantData svData, boolean linkOnStart)
    {
        mSvData = svData;
        mSvLinkOnStart = linkOnStart;
    }

    public final StructuralVariantData getStructuralVariantData() { return mSvData; }
    public boolean svLinkOnStart() { return mSvLinkOnStart; }

    public boolean matchesSegment(SegmentSupport segment, boolean isStart)
    {
        return isStart ? SegStart.equals(segment.toString()) : SegEnd.equals(segment.toString());
    }

    public static boolean isSvSegment(final SegmentSupport segment)
    {
        return (segment == SegmentSupport.BND || segment == SegmentSupport.INV || segment == SegmentSupport.INS
                || segment == SegmentSupport.DEL || segment == SegmentSupport.DUP || segment == SegmentSupport.SGL
                || segment == SegmentSupport.MULTIPLE);
     }

    public boolean matchesSV(boolean isStart)
    {
        return isStart ? isSvSegment(SegmentSupport.valueOf(SegStart)) : isSvSegment(SegmentSupport.valueOf(SegEnd));
    }

    public final String asString() { return String.format("id=%s pos=%s:%d", mId, Chromosome, StartPos); }

}
