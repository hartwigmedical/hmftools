package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_EVIDENCE_MAP_QUAL;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.QualityConfig;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.read.ExpandedBasesFactory;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.read.NumberEvents;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.VariantTier;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements VariantHotspot
{
    public final String Sample;
    public final VariantTier Tier;
    public final boolean Realign;
    public final int MaxCoverage;

    private final VariantHotspot mVariant;
    private final ReadContext mReadContext;
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;
    private final boolean mIsIndel;

    private int mFull;
    private int mPartial;
    private int mCore;
    private int mAlt;
    private int mRealigned;
    private int mReference;
    private int mCoverage;

    private int mLengthened;
    private int mShortened;

    private int mFullQuality;
    private int mPartialQuality;
    private int mCoreQuality;
    private int mAltQuality;
    private int mRealignedQuality;
    private int mReferenceQuality;
    private int mTotalQuality;

    private int mForwardStrand;
    private int mReverseStrand;

    private double mJitterPenalty;

    private int mImproperPair;

    private int mRawDepth;
    private int mRawAltSupport;
    private int mRawRefSupport;
    private int mRawAltBaseQuality;
    private int mRawRefBaseQuality;

    private int mLocalPhaseSet;

    public ReadContextCounter(
            final String sample, final VariantHotspot variant, final ReadContext readContext, final VariantTier tier,
            final int maxCoverage, final int minNumberOfEvents, boolean realign)
    {
        Sample = sample;
        Tier = tier;
        Realign = realign;
        MaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;

        mReadContext = readContext;
        mVariant = variant;
        mIsMnv = variant.isMNV();
        mIsIndel = variant.isIndel();

        mFull = 0;
        mPartial = 0;
        mCore = 0;
        mAlt = 0;
        mRealigned = 0;
        mReference = 0;
        mCoverage = 0;

        mLengthened = 0;
        mShortened = 0;

        mFullQuality = 0;
        mPartialQuality = 0;
        mCoreQuality = 0;
        mAltQuality = 0;
        mRealignedQuality = 0;
        mReferenceQuality = 0;
        mTotalQuality = 0;

        mForwardStrand = 0;
        mReverseStrand = 0;

        mJitterPenalty = 0;

        mImproperPair = 0;

        mRawDepth = 0;
        mRawAltSupport = 0;
        mRawRefSupport = 0;
        mRawAltBaseQuality = 0;
        mRawRefBaseQuality = 0;

        mLocalPhaseSet = 0;
    }

    public VariantHotspot variant() { return mVariant; }
    public ReadContext readContext() { return mReadContext; }

    @Override
    public String chromosome() { return mVariant.chromosome(); }

    @Override
    public int position() { return mVariant.position(); }

    @Override
    public String ref() { return mVariant.ref(); }

    @Override
    public String alt() { return mVariant.alt(); }

    public int altSupport() { return mFull + mPartial + mCore + mAlt + mRealigned; }

    public int refSupport() { return mReference; }
    public int coverage() { return mCoverage; }
    public int depth() { return mCoverage; }

    public double vaf() { return af(altSupport()); }

    private double af(double support)
    {
        return mCoverage == 0 ? 0d : support / mCoverage;
    }

    public int tumorQuality()
    {
        int tumorQuality = mFullQuality + mPartialQuality;
        return Math.max(0, tumorQuality - (int) mJitterPenalty);
    }

    public int[] counts()
    {
        return new int[] { mFull, mPartial, mCore, mRealigned, mAlt, mReference, mCoverage };
    }

    public int[] jitter()
    {
        return new int[] { mShortened, mLengthened, qualityJitterPenalty() };
    }

    public double strandBias()
    {
        double total = mForwardStrand + mReverseStrand;
        return total > 0 ? mForwardStrand / total : 0;
    }

    public int[] quality()
    {
        return new int[] { mFullQuality, mPartialQuality, mCoreQuality, mRealignedQuality, mAltQuality, mReferenceQuality, mTotalQuality };
    }

    public int improperPair() { return mImproperPair; }

    public int rawDepth() { return mRawDepth; }
    public int rawAltSupport() { return mRawAltSupport; }
    public int rawRefSupport() { return mRawRefSupport; }
    public int rawAltBaseQuality() { return mRawAltBaseQuality; }
    public int rawRefBaseQuality() { return mRawRefBaseQuality; }

    public void setLocalPhaseSet(int lps) { mLocalPhaseSet = lps; }
    public int localPhaseSet() { return mLocalPhaseSet; }

    public boolean exceedsMaxCoverage() { return mCoverage >= MaxCoverage; }

    public String toString()
    {
        return String.format("var(%s:%d %s>%s) core(%s) counts(f=%d p=%d c=%d)",
                mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(),
                mReadContext.toString(), mFull, mPartial, mCore);
    }

    public boolean processRead(final SAMRecord record, final SageConfig sageConfig, final QualityCalculator qualityCalc, final int rawNumberOfEvents)
    {
        if(exceedsMaxCoverage())
            return false;

        if(!Tier.equals(VariantTier.HOTSPOT) && record.getMappingQuality() < DEFAULT_EVIDENCE_MAP_QUAL)
            return false;

        final RawContext rawContext = RawContext.create(mVariant, record);

        if(rawContext.ReadIndexInSkipped)
            return false;

        int readIndex = rawContext.ReadIndex;
        boolean baseDeleted = rawContext.ReadIndexInDelete;

        mRawDepth += rawContext.DepthSupport ? 1 : 0;
        mRawAltSupport += rawContext.AltSupport ? 1 : 0;
        mRawRefSupport += rawContext.RefSupport ? 1 : 0;
        mRawAltBaseQuality += rawContext.AltQuality;
        mRawRefBaseQuality += rawContext.RefQuality;

        if(readIndex < 0)
            return false;

        boolean covered = mReadContext.isCoreCovered(readIndex, record.getReadBases());

        if(!covered)
            return false;

        final QualityConfig qualityConfig = sageConfig.Quality;

        int numberOfEvents = mIsMnv ?
                max(mMinNumberOfEvents, NumberEvents.numberOfEventsWithMNV(rawNumberOfEvents, mVariant.ref(), mVariant.alt()))
                : max(mMinNumberOfEvents, rawNumberOfEvents);

        double quality = qualityCalc.calculateQualityScore(this, readIndex, record, numberOfEvents);

        // Check if FULL, PARTIAL, OR CORE
        if(!baseDeleted)
        {
            boolean wildcardMatchInCore = mVariant.isSNV() && mReadContext.microhomology().isEmpty();

            IndexedBases readBases = record.getCigar().containsOperator(CigarOperator.N) ?
                    ExpandedBasesFactory.expand(position(), readIndex, record) :
                    new IndexedBases(position(), readIndex, record.getReadBases());

            int nonIndelLength = mIsIndel ? 0 : alt().length();

            final ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                    readBases, wildcardMatchInCore, record.getBaseQualities(), nonIndelLength);

            if(!match.equals(ReadContextMatch.NONE))
            {
                switch(match)
                {
                    case FULL:
                        mFull++;
                        mFullQuality += quality;
                        break;

                    case PARTIAL:
                        mPartial++;
                        mPartialQuality += quality;
                        break;

                    case CORE:
                        mCore++;
                        mCoreQuality += quality;
                        break;
                }

                mCoverage++;
                mTotalQuality += quality;

                countStrandedness(record);
                checkImproperCount(record);
                return true;
            }
        }

        // Check if REALIGNED
        final RealignedContext realignment = realignmentContext(Realign, readIndex, record);
        final RealignedType realignmentType = realignment.Type;
        if(realignmentType.equals(RealignedType.EXACT))
        {
            mRealigned++;
            mRealignedQuality += quality;
            mCoverage++;
            mTotalQuality += quality;
            return true;
        }

        if(realignmentType.equals(RealignedType.NONE) && rawContext.ReadIndexInSoftClip)
        {
            return false;
        }

        mCoverage++;
        mTotalQuality += quality;
        if(rawContext.RefSupport)
        {
            mReference++;
            mReferenceQuality += quality;
        }
        else if(rawContext.AltSupport)
        {
            mAlt++;
            mAltQuality++;
            countStrandedness(record);
        }

        // Jitter Penalty
        switch(realignmentType)
        {
            case LENGTHENED:
                mJitterPenalty += qualityConfig.jitterPenalty(realignment.RepeatCount);
                mLengthened++;
                break;
            case SHORTENED:
                mJitterPenalty += qualityConfig.jitterPenalty(realignment.RepeatCount);
                mShortened++;
                break;
        }

        return false;
    }

    private void countStrandedness(final SAMRecord record)
    {
        if(record.getFirstOfPairFlag())
            mForwardStrand++;
        else
            mReverseStrand++;
    }

    private RealignedContext realignmentContext(boolean realign, int readIndex, SAMRecord record)
    {
        if(!realign)
            return new RealignedContext(RealignedType.NONE, 0);

        int index = mReadContext.readBasesPositionIndex();
        int leftIndex = mReadContext.readBasesLeftCentreIndex();
        int rightIndex = mReadContext.readBasesRightCentreIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        int indelLength = indelLength(record);

        return Realigned.realignedAroundIndex(mReadContext,
                readIndex,
                record.getReadBases(),
                Math.max(indelLength + Math.max(leftOffset, rightOffset), Realigned.MAX_REPEAT_SIZE));
    }

    private int qualityJitterPenalty() { return (int) mJitterPenalty; }

    private void checkImproperCount(final SAMRecord record)
    {
        if(!record.getReadPairedFlag() || !record.getProperPairFlag())
        {
            mImproperPair++;
        }
    }

    private int indelLength(final SAMRecord record)
    {
        int result = 0;
        for(CigarElement cigarElement : record.getCigar())
        {
            switch(cigarElement.getOperator())
            {
                case I:
                case D:
                    result += cigarElement.getLength();
            }

        }

        return result;
    }
}
