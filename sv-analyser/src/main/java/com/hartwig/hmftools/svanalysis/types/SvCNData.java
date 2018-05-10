package com.hartwig.hmftools.svanalysis.types;

public class SvCNData {

    final private int mId;
    final private String mChromosome;
    final private long mStartPos;
    final private long mEndPos;
    final private double mCopyNumber;
    final private String mSegStart;
    final private String mSegEnd;
    final private String mMethod;
    final private double mBaf;

    public SvCNData(final int Id, final String Chromosome, final long StartPos, final long EndPos,
            final double CopyNumber, final String SegStart, final String SegEnd, final String Method, final double Baf)
    {
        mId = Id;
        mChromosome = Chromosome;
        mStartPos = StartPos;
        mEndPos = EndPos;
        mCopyNumber = CopyNumber;
        mSegStart = SegStart;
        mSegEnd = SegEnd;
        mMethod = Method;
        mBaf = Baf;
    }

    public int getId() { return mId; }
    public final String getChromosome() { return mChromosome; }
    public long getStartPos() { return mStartPos; }
    public long getEndPos() { return mEndPos; }
    public double getCopyNumber() { return mCopyNumber; }
    public final String getSegStart() { return mSegStart; }
    public final String getSegEnd() { return mSegEnd; }
    public final String getMethod() { return mMethod; }
    public double getBaf() { return mBaf; }

    public final String asString() { return String.format("id=%s pos=%s:%d", mId, mChromosome, mStartPos); }

}
