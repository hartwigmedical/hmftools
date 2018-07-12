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
    final private int mBafCount;
    final private double mObservedBaf;
    final private double mActualBaf;

    public static String CN_SEG_NONE = "NONE";
    public static String CN_SEG_UNKNOWN = "UNKNOWN";
    public static String CN_SEG_TELOMERE = "TELOMERE";
    public static String CN_SEG_CENTROMERE = "CENTROMERE";

    public SvCNData(final int Id, final String Chromosome, long StartPos, long EndPos,
            final String SegStart, final String SegEnd, int bafCount, double observedBaf, double actualBaf,
            double copyNumber, final String method)
    {
        mId = Id;
        mChromosome = Chromosome;
        mStartPos = StartPos;
        mEndPos = EndPos;
        mCopyNumber = copyNumber;
        mSegStart = SegStart;
        mSegEnd = SegEnd;
        mMethod = method;
        mBafCount = bafCount;
        mObservedBaf = observedBaf;
        mActualBaf = actualBaf;
    }

    public int id() { return mId; }
    public final String chromosome() { return mChromosome; }
    public long startPos() { return mStartPos; }
    public long endPos() { return mEndPos; }
    public double copyNumber() { return mCopyNumber; }
    public final String segStart() { return mSegStart; }
    public final String segEnd() { return mSegEnd; }
    public final String method() { return mMethod; }
    public int bafCount() { return mBafCount; }
    public double observedBaf() { return mObservedBaf; }
    public double actualBaf() { return mActualBaf; }

    public final String asString() { return String.format("id=%s pos=%s:%d", mId, mChromosome, mStartPos); }

}
