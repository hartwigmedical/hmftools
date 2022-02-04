package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_EVIDENCE_MAP_QUAL;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.UNRELATED;

import java.util.Comparator;
import java.util.List;

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

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements VariantHotspot
{
    public final VariantTier Tier;
    public final boolean Realign;
    public final int MaxCoverage;

    private final int mId;

    private final VariantHotspot mVariant;
    private final ReadContext mReadContext;
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;
    private final boolean mIsIndel;

    private final int[] mQualities;
    private final int[] mCounts;

    private int mLengthened;
    private int mShortened;

    private int mForwardStrand;
    private int mReverseStrand;

    private double mJitterPenalty;

    private int mImproperPair;

    private int mRawDepth;
    private int mRawAltSupport;
    private int mRawRefSupport;
    private int mRawAltBaseQuality;
    private int mRawRefBaseQuality;

    private List<Integer> mLocalPhaseSets;
    private List<double[]> mLpsCounts;

    public static final int RC_FULL = 0;
    public static final int RC_PARTIAL = 1;
    public static final int RC_CORE = 2;
    public static final int RC_REALIGNED = 3;
    public static final int RC_ALT = 4;
    public static final int RC_REF = 5;
    public static final int RC_TOTAL = 6;
    public static final int RC_MAX = RC_TOTAL + 1;

    public ReadContextCounter(
            final int id, final VariantHotspot variant, final ReadContext readContext, final VariantTier tier,
            final int maxCoverage, final int minNumberOfEvents, boolean realign)
    {
        mId = id;

        Tier = tier;
        Realign = realign;
        MaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;

        mReadContext = readContext;
        mVariant = variant;
        mIsMnv = variant.isMNV();
        mIsIndel = variant.isIndel();

        mQualities = new int[RC_MAX];
        mCounts = new int[RC_MAX];

        mLengthened = 0;
        mShortened = 0;

        mForwardStrand = 0;
        mReverseStrand = 0;

        mJitterPenalty = 0;

        mImproperPair = 0;

        mRawDepth = 0;
        mRawAltSupport = 0;
        mRawRefSupport = 0;
        mRawAltBaseQuality = 0;
        mRawRefBaseQuality = 0;

        mLocalPhaseSets = null;
        mLpsCounts = null;
    }

    public int id() { return mId; }
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

    public int altSupport() { return mCounts[RC_FULL] + mCounts[RC_PARTIAL] + mCounts[RC_CORE] + mCounts[RC_ALT] + mCounts[RC_REALIGNED]; }

    public int refSupport() { return mCounts[RC_REF]; }
    public int coverage() { return mCounts[RC_TOTAL]; }
    public int depth() { return mCounts[RC_TOTAL]; }

    public double vaf() { return af(altSupport()); }

    private double af(double support)
    {
        return mCounts[RC_TOTAL] == 0 ? 0d : support / mCounts[RC_TOTAL];
    }

    public int tumorQuality()
    {
        int tumorQuality = mQualities[RC_FULL] + mQualities[RC_PARTIAL];
        return Math.max(0, tumorQuality - (int) mJitterPenalty);
    }

    public int[] counts() { return mCounts; }
    public int[] quality() { return mQualities; }

    public int[] jitter()
    {
        return new int[] { mShortened, mLengthened, qualityJitterPenalty() };
    }

    public double strandBias()
    {
        double total = mForwardStrand + mReverseStrand;
        return total > 0 ? mForwardStrand / total : 0;
    }

    public int improperPair() { return mImproperPair; }

    public int rawDepth() { return mRawDepth; }
    public int rawAltSupport() { return mRawAltSupport; }
    public int rawRefSupport() { return mRawRefSupport; }
    public int rawAltBaseQuality() { return mRawAltBaseQuality; }
    public int rawRefBaseQuality() { return mRawRefBaseQuality; }

    public void addLocalPhaseSet(int lps, int readCount, double allocCount)
    {
        if(mLocalPhaseSets == null)
        {
            mLocalPhaseSets = Lists.newArrayList();
            mLpsCounts = Lists.newArrayList();
        }

        mLocalPhaseSets.add(lps);
        mLpsCounts.add(new double[] { readCount, allocCount } );
    }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<double[]> lpsCounts() { return mLpsCounts; }

    public boolean exceedsMaxCoverage() { return mCounts[RC_TOTAL] >= MaxCoverage; }

    public String toString()
    {
        return String.format("id(%d) var(%s:%d %s>%s) core(%s) counts(f=%d p=%d c=%d)",
                mId, mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt(),
                mReadContext.toString(), mCounts[RC_FULL], mCounts[RC_PARTIAL], mCounts[RC_CORE]);
    }

    public ReadMatchType processRead(final SAMRecord record, final SageConfig sageConfig, final QualityCalculator qualityCalc, final int rawNumberOfEvents)
    {
        if(exceedsMaxCoverage())
            return UNRELATED;

        if(!Tier.equals(VariantTier.HOTSPOT) && record.getMappingQuality() < DEFAULT_EVIDENCE_MAP_QUAL)
            return UNRELATED;

        final RawContext rawContext = RawContext.create(mVariant, record);

        if(rawContext.ReadIndexInSkipped)
            return UNRELATED;

        int readIndex = rawContext.ReadIndex;
        boolean baseDeleted = rawContext.ReadIndexInDelete;

        mRawDepth += rawContext.DepthSupport ? 1 : 0;
        mRawAltSupport += rawContext.AltSupport ? 1 : 0;
        mRawRefSupport += rawContext.RefSupport ? 1 : 0;
        mRawAltBaseQuality += rawContext.AltQuality;
        mRawRefBaseQuality += rawContext.RefQuality;

        if(readIndex < 0)
            return UNRELATED;

        boolean covered = mReadContext.isCoreCovered(readIndex, record.getReadBases());

        if(!covered)
            return UNRELATED;

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
                        mCounts[RC_FULL]++;
                        mQualities[RC_FULL] += quality;
                        break;

                    case PARTIAL:
                        mCounts[RC_PARTIAL]++;
                        mQualities[RC_PARTIAL] += quality;
                        break;

                    case CORE:
                        ++mCounts[RC_CORE];
                        mQualities[RC_CORE] += quality;
                        break;
                }

                ++mCounts[RC_TOTAL];
                mQualities[RC_TOTAL] += quality;

                countStrandedness(record);
                checkImproperCount(record);
                return SUPPORT;
            }
        }

        // Check if REALIGNED
        final RealignedContext realignment = realignmentContext(Realign, readIndex, record);
        final RealignedType realignmentType = realignment.Type;
        if(realignmentType.equals(RealignedType.EXACT))
        {
            mCounts[RC_REALIGNED]++;
            mQualities[RC_REALIGNED] += quality;

            mCounts[RC_TOTAL]++;
            mQualities[RC_TOTAL] += quality;
            return SUPPORT;
        }

        if(realignmentType.equals(RealignedType.NONE) && rawContext.ReadIndexInSoftClip)
        {
            return UNRELATED;
        }

        ReadMatchType matchType = UNRELATED;

        mCounts[RC_TOTAL]++;
        mQualities[RC_TOTAL] += quality;

        if(rawContext.RefSupport)
        {
            mCounts[RC_REF]++;
            mQualities[RC_REF] += quality;
            matchType = NO_SUPPORT;
        }
        else if(rawContext.AltSupport)
        {
            mCounts[RC_ALT]++;
            mQualities[RC_ALT] += quality;
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

        return matchType;
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
