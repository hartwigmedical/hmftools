package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.purple.PurpleConstants.MAX_DIPLOID_COPY_NUMBER;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_DIPLOID_COPY_NUMBER;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class ObservedRegion implements GenomeRegion
{
    private final String mChromosome;
    private int mPosStart;
    private int mPosEnd;

    private boolean mRatioSupport;

    private SegmentSupport mSupport;

    private int mBafCount;
    private double mObservedBAF;
    private int mDepthWindowCount;
    private double mObservedTumorRatio;
    private double mObservedNormalRatio;
    private double mUnnormalisedObservedNormalRatio;

    private GermlineStatus mGermlineStatus;
    private boolean mSvCluster;
    private double mGcContent;

    private int mMinStart;
    private int mMaxStart;

    // fitted region fields
    private double mMinorAlleleCopyNumberDeviation;
    private double mMajorAlleleCopyNumberDeviation;
    private double mDeviationPenalty;
    private double mEventPenalty;
    private double mRefNormalisedCopyNumber;
    private double mTumorCopyNumber;
    private double mTumorBAF;
    private double mFittedTumorCopyNumber;
    private double mFittedBAF;

    public ObservedRegion(final String chromosome, final int posStart, final int posEnd, final boolean ratioSupport,
            final SegmentSupport support, final int bafCount, final double observedBAF, final int depthWindowCount,
            final double observedTumorRatio, final double observedNormalRatio, final double unnormalisedObservedNormalRatio,
            final GermlineStatus germlineStatus, final boolean svCluster, final double gcContent, final int minStart, final int maxStart,
            double minorAlleleCopyNumberDeviation, double majorAlleleCopyNumberDeviation, double deviationPenalty, double eventPenalty,
            double refNormalisedCopyNumber, double tumorCopyNumber, double tumorBAF, double fittedTumorCopyNumber, double fittedBAF)
    {
        mChromosome = chromosome;
        mPosStart = posStart;
        mPosEnd = posEnd;

        mRatioSupport = ratioSupport;
        mSupport = support;
        mBafCount = bafCount;
        mObservedBAF = observedBAF;
        mDepthWindowCount = depthWindowCount;
        mObservedTumorRatio = observedTumorRatio;
        mObservedNormalRatio = observedNormalRatio;
        mUnnormalisedObservedNormalRatio = unnormalisedObservedNormalRatio;
        mGermlineStatus = germlineStatus;
        mSvCluster = svCluster;
        mGcContent = gcContent;
        mMinStart = minStart;
        mMaxStart = maxStart;

        mMinorAlleleCopyNumberDeviation = minorAlleleCopyNumberDeviation;
        mMajorAlleleCopyNumberDeviation = majorAlleleCopyNumberDeviation;
        mDeviationPenalty = deviationPenalty;
        mEventPenalty = eventPenalty;
        mRefNormalisedCopyNumber = refNormalisedCopyNumber;
        mTumorCopyNumber = tumorCopyNumber;
        mTumorBAF = tumorBAF;
        mFittedTumorCopyNumber = fittedTumorCopyNumber;
        mFittedBAF = fittedBAF;
    }

    public static ObservedRegion from(final ObservedRegion other)
    {
        return new ObservedRegion(
                other.chromosome(), other.start(), other.end(),
                other.ratioSupport(), other.support(), other.bafCount(), other.observedBAF(), other.depthWindowCount(),
                other.observedTumorRatio(), other.observedNormalRatio(), other.unnormalisedObservedNormalRatio(), other.germlineStatus(),
                other.svCluster(), other.gcContent(), other.minStart(), other.maxStart(), other.minorAlleleCopyNumberDeviation(),
                other.majorAlleleCopyNumberDeviation(), other.deviationPenalty(), other.eventPenalty(),
                other.refNormalisedCopyNumber(), other.tumorCopyNumber(), other.tumorBAF(), other.fittedTumorCopyNumber(), other.fittedBAF());
    }

    public ObservedRegion(
            final String chromosome, final int posStart, final int posEnd, final boolean ratioSupport,
            final SegmentSupport support, final int bafCount, final double observedBAF, final int depthWindowCount,
            final double observedTumorRatio, final double observedNormalRatio, final double unnormalisedObservedNormalRatio,
            final GermlineStatus germlineStatus, final boolean svCluster, final double gcContent, final int minStart, final int maxStart)
    {
        this(chromosome, posStart, posEnd, ratioSupport, support, bafCount, observedBAF, depthWindowCount, observedTumorRatio,
                observedNormalRatio, unnormalisedObservedNormalRatio, germlineStatus, svCluster, gcContent, minStart, maxStart,
                0, 0, 0, 0, 0,
                0, 0,0, 0);
    }
    @Override
    public String chromosome() { return mChromosome; }

    @Override
    public int start() { return mPosStart; }
    public int end() { return mPosEnd; }

    public void setStart(int start) { mPosStart = start; }
    public void setEnd(int end) { mPosEnd = end; }

    public int minStart() { return mMinStart; }
    public int maxStart() { return mMaxStart; }

    public void setMinStart(int start) { mMinStart = start; }
    public void setMaxStart(int start) { mMaxStart = start; }

    public boolean ratioSupport() { return mRatioSupport; }
    public void setRatioSupport(boolean ratioSupport) { mRatioSupport = ratioSupport; }

    public SegmentSupport support() { return mSupport; }
    public void setSupport(final SegmentSupport support) { mSupport = support; }

    public int bafCount() { return mBafCount; }
    public void setBafCount(int bafCount) { mBafCount = bafCount; }

    public double observedBAF() { return mObservedBAF; }
    public void setObservedBAF(double baf) { mObservedBAF = baf; }

    public int depthWindowCount() { return mDepthWindowCount; }
    public void setDepthWindowCount(int count) { mDepthWindowCount = count; }

    public double observedTumorRatio() { return mObservedTumorRatio; }
    public void setObservedTumorRatio(double ratio) { mObservedTumorRatio = ratio; }

    public double observedNormalRatio() { return mObservedNormalRatio; }
    public void setObservedNormalRatio(double ratio) { mObservedNormalRatio = ratio; }

    public double unnormalisedObservedNormalRatio() { return mUnnormalisedObservedNormalRatio; }

    public GermlineStatus germlineStatus() { return mGermlineStatus; }
    public void setGermlineStatus(final GermlineStatus status)
    {
        mGermlineStatus = status;
    }

    public boolean svCluster() { return mSvCluster; }

    public double gcContent() { return mGcContent; }
    public void setGcContent(double content) { mGcContent = content; }

    // fitted data and methods
    public double minorAlleleCopyNumberDeviation() { return mMinorAlleleCopyNumberDeviation; }
    public void setMinorAlleleCopyNumberDeviation(double deviation) { mMinorAlleleCopyNumberDeviation = deviation; }

    public double majorAlleleCopyNumberDeviation() { return mMajorAlleleCopyNumberDeviation; }
    public void setMajorAlleleCopyNumberDeviation(double deviation) { mMajorAlleleCopyNumberDeviation = deviation; }

    public double deviationPenalty() { return mDeviationPenalty; }
    public void setDeviationPenalty(double penalty) { mDeviationPenalty = penalty; }

    public double eventPenalty() { return mEventPenalty; }
    public void setEventPenalty(double penalty) { mEventPenalty = penalty; }

    public double refNormalisedCopyNumber() { return mRefNormalisedCopyNumber; }
    public void setRefNormalisedCopyNumber(double cn) { mRefNormalisedCopyNumber = cn; }

    public double tumorCopyNumber() { return mTumorCopyNumber; }
    public void setTumorCopyNumber(double cn) { mTumorCopyNumber = cn; }

    public double tumorBAF() { return mTumorBAF; }
    public void setTumorBAF(double baf) { mTumorBAF = baf; }

    public double fittedTumorCopyNumber() { return mFittedTumorCopyNumber; }
    public void setFittedTumorCopyNumber(double cn) { mFittedTumorCopyNumber = cn; }

    public double fittedBAF() { return mFittedBAF; }
    public void setFittedBAF(double baf) { mFittedBAF = baf; }

    public double minorAlleleCopyNumber() { return mTumorCopyNumber - majorAlleleCopyNumber(); }
    public double majorAlleleCopyNumber() { return mTumorBAF * mTumorCopyNumber; }

    public boolean isDiploid()
    {
        return Doubles.greaterOrEqual(majorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER)
                && Doubles.lessOrEqual(majorAlleleCopyNumber(), MAX_DIPLOID_COPY_NUMBER)
                && Doubles.greaterOrEqual(minorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER)
                && Doubles.lessOrEqual(minorAlleleCopyNumber(), MAX_DIPLOID_COPY_NUMBER);
    }

    public String toString()
    {
        return String.format("loc(%s:%d-%d) ratio(%s) seg(%s) status(%s) baf(%d obs=%.2f) cn(%.2f ref=%.2f baf=%.2f) fit(%.2f baf=%.2f) start(min=%d max=%d)",
                mChromosome, mPosStart, mPosEnd, mRatioSupport, mSupport, mGermlineStatus, mBafCount, mObservedBAF,
                mTumorCopyNumber, mRefNormalisedCopyNumber, mTumorBAF, mFittedTumorCopyNumber, mFittedBAF,
                mMinStart, mMaxStart);
    }
}
