package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.CORE;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.FULL;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.OTHER_ALT;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.PARTIAL;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REALIGNED;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REF;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_BASE_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.EVIDENCE_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MQ_RATIO_SMOOTHING;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.candidate.RefContextConsumer.ignoreSoftClipAdapter;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.ReadEdgeDistance.calcAdjustedVariantPosition;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.ALT_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.CHIMERIC;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.IN_SPLIT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAP_QUAL;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAX_COVERAGE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NON_CORE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.REF_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SOFT_CLIP;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.UNRELATED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.CORE_PARTIAL;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;
import static com.hartwig.hmftools.sage.filter.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.filter.FragmentCoords;
import com.hartwig.hmftools.sage.filter.StrandBiasData;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.read.NumberEvents;
import com.hartwig.hmftools.sage.read.SplitReadUtils;
import com.hartwig.hmftools.sage.vis.VariantVis;
import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter//  extends SimpleVariant
{
    private final int mId;
    private final VariantTier mTier;
    private final SimpleVariant mVariant;
    private final ReadContext mReadContext;
    private final SageConfig mConfig;
    private final QualityCalculator mQualityCalculator;
    private final String mSample;
    private final int mMaxCoverage;
    private final VariantVis mVariantVis;

    // local variant-related state
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;
    private final boolean mIsSnv;
    private final boolean mIsIndel;
    private final boolean mAllowWildcardMatchInCore;
    private final int mMaxCoreMismatches;
    private final BqrQualCache mBqrQualCache;

    // counts of various
    private final ReadSupportCounts mQualities;
    private final ReadSupportCounts mCounts;

    private final StrandBiasData mAltFragmentStrandBias;
    private final StrandBiasData mRefFragmentStrandBias;
    private final StrandBiasData mAltReadStrandBias;
    private final StrandBiasData mRefReadStrandBias;

    private final JitterData mJitterData;
    private int mImproperPairCount;

    private int mRawDepth;
    private int mRawAltSupport;
    private int mRawRefSupport;
    private int mRawAltBaseQualityTotal;
    private int mRawRefBaseQualityTotal;
    private int mRecalibratedBaseQualityTotal;
    private double mSupportAltBaseQualityTotal;

    private long mMapQualityTotal;
    private long mAltMapQualityTotal;
    private int mNmCountTotal;
    private int mAltNmCountTotal;

    private int mSoftClipInsertSupport;
    private int mMaxCandidateDeleteLength;
    private final ReadEdgeDistance mReadEdgeDistance;

    private List<Integer> mLocalPhaseSets;
    private List<Integer> mLpsCounts;
    private int[] mUmiTypeCounts;
    private FragmentLengthData mFragmentLengthData;
    private FragmentCoords mFragmentCoords;

    public ReadContextCounter(
            final int id, final SimpleVariant variant, final ReadContext readContext, final VariantTier tier,
            final int maxCoverage, final int minNumberOfEvents, final SageConfig config, final QualityCalculator qualityCalculator,
            final String sampleId)
    {
        mId = id;

        mTier = tier;
        mMaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;
        mSample = sampleId;
        mQualityCalculator = qualityCalculator;
        mConfig = config;

        mReadContext = readContext;
        mVariant = variant;

        mVariantVis = config.Visualiser.Enabled && config.Visualiser.processVariant(variant)
                ? new VariantVis(mConfig, mSample, mVariant, mReadContext, mTier) : null;

        // set local state to avoid testing on each read
        mIsMnv = variant.isMNV();
        mIsSnv = variant.isSNV();
        mIsIndel = variant.isIndel();
        mBqrQualCache = !mIsIndel ? new BqrQualCache(variant.position(), variant.alt()) : null;

        mAllowWildcardMatchInCore = mVariant.isSNV() && mReadContext.microhomology().isEmpty();

        mMaxCoreMismatches = mVariant.isIndel() && mVariant.alt().length() >= CORE_LOW_QUAL_MISMATCH_BASE_LENGTH ?
                mVariant.alt().length() / CORE_LOW_QUAL_MISMATCH_BASE_LENGTH : 0;

        mQualities = new ReadSupportCounts();
        mCounts = new ReadSupportCounts();

        mJitterData = new JitterData();

        mAltFragmentStrandBias = new StrandBiasData(true);
        mRefFragmentStrandBias = new StrandBiasData(false);
        mAltReadStrandBias = new StrandBiasData(true);
        mRefReadStrandBias = new StrandBiasData(false);

        mImproperPairCount = 0;

        mRawDepth = 0;
        mRawAltSupport = 0;
        mRawRefSupport = 0;
        mRawAltBaseQualityTotal = 0;
        mRawRefBaseQualityTotal = 0;
        mRecalibratedBaseQualityTotal = 0;
        mSupportAltBaseQualityTotal = 0;
        mSoftClipInsertSupport = 0;
        mMaxCandidateDeleteLength = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mNmCountTotal = 0;
        mAltNmCountTotal = 0;

        mReadEdgeDistance = new ReadEdgeDistance(calcAdjustedVariantPosition(mVariant.position(), indelLength()));

        mLocalPhaseSets = null;
        mLpsCounts = null;
        mUmiTypeCounts = null;
        mFragmentLengthData = mConfig.WriteFragmentLengths ? new FragmentLengthData() : null;
        mFragmentCoords = mConfig.Quality.HighBaseMode ? new FragmentCoords(REQUIRED_UNIQUE_FRAG_COORDS) : null;
    }

    public int id() { return mId; }
    public SimpleVariant variant() { return mVariant; }
    public ReadContext readContext() { return mReadContext; }
    public VariantTier tier() { return mTier; }
    public int indelLength() { return mVariant.isIndel() ? max(mVariant.alt().length(), mVariant.ref().length()) : 0; }
    public boolean isSnv() { return mIsSnv; }
    public boolean isIndel() { return mIsIndel; }
    public final BqrQualCache bqrQualCache() { return mBqrQualCache; }
    public String chromosome() { return mVariant.chromosome(); }
    public int position() { return mVariant.position(); }
    public String ref() { return mVariant.ref(); }
    public String alt() { return mVariant.alt(); }

    public int altSupport() { return mCounts.altSupport(); }
    public int strongAltSupport() { return mCounts.strongSupport(); }
    public int refSupport() { return mCounts.Ref; }
    public int depth() { return mCounts.Total; }

    public double vaf()
    {
        return mCounts.Total == 0 ? 0d : mCounts.altSupport() / (double)mCounts.Total;
    }

    public int tumorQuality()
    {
        int tumorQuality = mQualities.Full + mQualities.Partial + mQualities.Realigned;
        return Math.max(0, tumorQuality - mJitterData.penalty());
    }

    public int[] counts() { return mCounts.toArray(); }
    public int[] quality() { return mQualities.toArray(); }

    public int[] jitter()
    {
        return mJitterData.summary();
    }

    public StrandBiasData fragmentStrandBiasAlt() { return mAltFragmentStrandBias; }
    public StrandBiasData fragmentStrandBiasRef() { return mRefFragmentStrandBias; }
    public StrandBiasData readStrandBiasAlt() { return mAltReadStrandBias; }
    public StrandBiasData readStrandBiasRef() { return mRefReadStrandBias; }
    public int improperPairCount() { return mImproperPairCount; }

    public int rawDepth() { return mRawDepth; }
    public int rawAltSupport() { return mRawAltSupport; }
    public int rawRefSupport() { return mRawRefSupport; }
    public int rawAltBaseQualityTotal() { return mRawAltBaseQualityTotal; }
    public int rawRefBaseQualityTotal() { return mRawRefBaseQualityTotal; }
    public int recalibratedBaseQualityTotal() { return mRecalibratedBaseQualityTotal; }

    public long mapQualityTotal() { return mMapQualityTotal; }
    public long altMapQualityTotal() { return mAltMapQualityTotal; }

    public long nmCountTotal() { return mNmCountTotal; }
    public long altNmCountTotal() { return mAltNmCountTotal; }
    public ReadEdgeDistance readEdgeDistance() { return mReadEdgeDistance; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    @Nullable
    public VariantVis variantVis() { return mVariantVis; }

    public double averageAltBaseQuality()
    {
        // excludes realigned
        int supportCount = mCounts.Full + mCounts.Partial + mCounts.Core;
        return supportCount > 0 ? mSupportAltBaseQualityTotal / (double)supportCount : 0;
    }

    public int softClipInsertSupport() { return mSoftClipInsertSupport; }

    public void setMaxCandidateDeleteLength(int length) { mMaxCandidateDeleteLength = length; }
    public int maxCandidateDeleteLength() { return mMaxCandidateDeleteLength; }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<Integer> lpsCounts() { return mLpsCounts; }

    public int[] umiTypeCounts() { return mUmiTypeCounts; }
    public FragmentLengthData fragmentLengths() { return mFragmentLengthData; }

    public boolean exceedsMaxCoverage() { return mCounts.Total >= mMaxCoverage; }

    public boolean belowMinFragmentCoords() { return mFragmentCoords != null && !mFragmentCoords.atCapacity(); }

    public String toString()
    {
        return format("id(%d) var(%s) core(%s) counts(f=%d p=%d c=%d)",
                mId, varString(), mReadContext.toString(), mCounts.Full, mCounts.Partial, mCounts.Core);
    }

    public String varString()
    {
        return format("%s:%d %s>%s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
    }

    public enum MatchType implements Comparable<MatchType>
    {
        FULL(0),
        PARTIAL(1),
        CORE(2),
        REALIGNED(3),
        CORE_PARTIAL(4),
        ALT(5),
        REF(6),
        NONE(7);

        public final int SortKey;

        MatchType(int sortKey)
        {
            SortKey = sortKey;
        }
    }

    public ReadMatchType processRead(final SAMRecord record, int numberOfEvents, @Nullable final FragmentData fragmentData)
    {
        if(exceedsMaxCoverage())
            return MAX_COVERAGE;

        if(mTier != VariantTier.HOTSPOT && record.getMappingQuality() < EVIDENCE_MIN_MAP_QUAL)
        {
            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return MAP_QUAL;
        }

        if(mConfig.Quality.HighBaseMode && isChimericRead(record))
        {
            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return CHIMERIC;
        }

        RawContext rawContext = RawContext.create(mVariant, record);

        if(mConfig.Quality.HighBaseMode && rawContext.ReadIndexInSoftClip)
        {
            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return SOFT_CLIP;
        }

        if(rawContext.ReadIndex < 0 && !ignoreSoftClipAdapter(record))
        {
            if(rawContext.DepthSupport || rawContext.AltSupport || rawContext.RefSupport)
            {
                SG_LOGGER.error("rawContext missing readIndex but with support(depth={} ref={} alt={}",
                        rawContext.DepthSupport, rawContext.AltSupport, rawContext.RefSupport);
            }

            // search for a core match within soft-clipped bases, checking if a proixmate DEL may explain the soft-clipping
            rawContext = createRawContextFromCoreMatch(record);

            if(rawContext.ReadIndex < 0)
            {
                addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
                return UNRELATED;
            }
        }

        if(rawContext.ReadIndexInSkipped)
        {
            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return IN_SPLIT;
        }

        int readIndex = rawContext.ReadIndex;
        boolean baseDeleted = rawContext.ReadIndexInDelete;

        mRawDepth += rawContext.DepthSupport ? 1 : 0;

        if(rawContext.ReadIndexInSoftClip && rawContext.AltSupport)
            ++mSoftClipInsertSupport;

        boolean covered = mReadContext.indexedBases().isCoreCovered(readIndex, record.getReadBases().length);

        if(!covered)
        {
            registerRawSupport(rawContext, 0);
            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return NON_CORE;
        }

        double adjustedNumOfEvents = numberOfEvents;

        if(mIsMnv)
            adjustedNumOfEvents = NumberEvents.calcWithMnvRaw(numberOfEvents, mVariant.ref(), mVariant.alt());

        if(max(mVariant.ref().length(), mVariant.alt().length()) <= SC_READ_EVENTS_FACTOR)
        {
            // penalise variants except long INDELs for their soft-clipped bases
            adjustedNumOfEvents += NumberEvents.calcSoftClipAdjustment(record);
        }

        adjustedNumOfEvents = max(mMinNumberOfEvents, adjustedNumOfEvents);

        double rawBaseQuality = mQualityCalculator.rawBaseQuality(this, readIndex, record);

        if(mConfig.Quality.HighBaseMode && rawBaseQuality < mConfig.Quality.HighBaseQualLimit)
        {
            if(rawContext.AltSupport)
                countAltSupportMetrics(record, fragmentData);

            addVariantVisRecord(record, MatchType.NONE, null, fragmentData);
            return UNRELATED;
        }

        QualityCalculator.QualityScores qualityScores = mQualityCalculator.calculateQualityScores(
                this, readIndex, record, adjustedNumOfEvents, rawBaseQuality);

        double quality = qualityScores.ModifiedQuality;

        MatchType matchType = MatchType.NONE;

        // Check if FULL, PARTIAL, OR CORE
        if(!baseDeleted)
        {
            final ReadContextMatch match = determineReadContextMatch(record, readIndex, true);

            if(match != NONE && match != ReadContextMatch.CORE_PARTIAL)
            {
                VariantReadSupport readSupport = null;

                switch(match)
                {
                    case FULL:
                        readSupport = FULL;
                        matchType = MatchType.FULL;
                        break;

                    case PARTIAL:
                        readSupport = PARTIAL;
                        matchType = MatchType.PARTIAL;
                        break;

                    case CORE:
                        readSupport = CORE;
                        matchType = MatchType.CORE;
                        break;
                }

                registerReadSupport(record, readSupport, quality);

                mMapQualityTotal += record.getMappingQuality();
                mAltMapQualityTotal += record.getMappingQuality();
                mNmCountTotal += numberOfEvents;
                mAltNmCountTotal += numberOfEvents;

                mSupportAltBaseQualityTotal += rawBaseQuality;

                registerRawSupport(rawContext, qualityScores.RecalibratedBaseQuality);

                if(rawContext.AltSupport)
                    mReadEdgeDistance.update(record, fragmentData, true);

                addVariantVisRecord(record, matchType, qualityScores, fragmentData);
                logReadEvidence(record, matchType, readIndex, quality);

                /*
                if(SG_LOGGER.isTraceEnabled() && sampleId != null)
                {
                    qualityCalc.logReadQualCalcs(this, readIndex, record, adjustedNumOfEvents);
                }
                */

                if(rawContext.AltSupport || readSupport != null)
                    countAltSupportMetrics(record, fragmentData);

                checkImproperCount(record);
                return ALT_SUPPORT;
            }
            else if(match == ReadContextMatch.CORE_PARTIAL)
            {
                // if the core is partly overlapped then back out any attribution to a ref match
                rawContext.updateSupport(false, rawContext.AltSupport);
            }
        }

        boolean canRealign = abs(mVariant.indelLength()) >= REALIGN_READ_MIN_INDEL_LENGTH || readHasIndelInCore(record);
        RealignedContext realignment = canRealign ? checkRealignment(record) : RealignedContext.NONE;

        if(realignment.Type == EXACT)
        {
            registerReadSupport(record, REALIGNED, quality);

            mMapQualityTotal += record.getMappingQuality();
            mAltMapQualityTotal += record.getMappingQuality();
            mNmCountTotal += numberOfEvents;
            mAltNmCountTotal += numberOfEvents;

            addVariantVisRecord(record, MatchType.REALIGNED, qualityScores, fragmentData);
            logReadEvidence(record, MatchType.REALIGNED, readIndex,quality);
            rawContext.updateSupport(false, rawContext.AltSupport);
            registerRawSupport(rawContext, qualityScores.RecalibratedBaseQuality);

            return ALT_SUPPORT;
        }
        else if(realignment.Type == CORE_PARTIAL)
        {
            matchType = MatchType.CORE_PARTIAL;
            rawContext.updateSupport(false, false);
        }

        registerRawSupport(rawContext, qualityScores.RecalibratedBaseQuality);

        // switch back to the old method to test for jitter
        RealignedContext jitterRealign = Realignment.realignedAroundIndex(mReadContext, readIndex, record.getReadBases(), getMaxRealignDistance(record));

        if(rawContext.ReadIndexInSoftClip && !rawContext.AltSupport)
        {
            if(jitterRealign.Type != LENGTHENED && jitterRealign.Type != SHORTENED)
            {
                addVariantVisRecord(record, MatchType.NONE, qualityScores, fragmentData);
                return SOFT_CLIP;
            }
        }

        ReadMatchType readMatchType = UNRELATED;

        mMapQualityTotal += record.getMappingQuality();
        mNmCountTotal += numberOfEvents;

        VariantReadSupport readSupport = null;

        if(rawContext.RefSupport)
        {
            readSupport = REF;
            readMatchType = REF_SUPPORT;

            mRefFragmentStrandBias.registerFragment(record);
            mRefReadStrandBias.registerRead(record, fragmentData, mVariant);

            mReadEdgeDistance.update(record, fragmentData, false);
        }
        else if(rawContext.AltSupport)
        {
            readSupport = OTHER_ALT;

            mAltMapQualityTotal += record.getMappingQuality();
            mAltNmCountTotal += numberOfEvents;

            countAltSupportMetrics(record, fragmentData);
        }

        registerReadSupport(record, readSupport, quality);

        // add to jitter penalty as a function of the number of repeats found
        mJitterData.update(jitterRealign, mConfig.Quality);

        if(rawContext.RefSupport)
            matchType = MatchType.REF;
        else if(rawContext.AltSupport)
            matchType = MatchType.ALT;

        addVariantVisRecord(record, matchType, qualityScores, fragmentData);
        logReadEvidence(record, matchType, readIndex, quality);

        return readMatchType;
    }

    private ReadContextMatch determineReadContextMatch(final SAMRecord record, int readIndex, boolean allowCoreVariation)
    {
        ReadIndexBases readIndexBases;

        if(record.getCigar().containsOperator(CigarOperator.N))
        {
            readIndexBases = SplitReadUtils.expandSplitRead(readIndex, record);
        }
        else
        {
            readIndexBases = new ReadIndexBases(readIndex, record.getReadBases());
        }

        final ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                readIndexBases, record.getBaseQualities(),
                allowCoreVariation ? mAllowWildcardMatchInCore : false,
                allowCoreVariation ? mMaxCoreMismatches : 0);

        return match;
    }

    private void registerReadSupport(final SAMRecord record, @Nullable final VariantReadSupport support, final double quality)
    {
        mCounts.addSupport(support, 1);
        mQualities.addSupport(support, (int)quality);

        boolean supportsVariant = support != null
                && (support == FULL || support == PARTIAL || support == CORE || support == REALIGNED);

        if(mConfig.TrackUMIs)
        {
            countUmiType(record, supportsVariant);
        }

        if(mFragmentLengthData != null && (support == REF || supportsVariant))
        {
            mFragmentLengthData.addLength(abs(record.getInferredInsertSize()), supportsVariant);
        }
    }

    private void countUmiType(final SAMRecord record, final boolean supportsVariant)
    {
        if(mUmiTypeCounts == null)
        {
            // 3 total depth values followed by the 3 variant support values
            mUmiTypeCounts = new int[UMI_TYPE_COUNT];
        }

        UmiReadType umiReadType = extractUmiType(record);

        // add to total and variant support if applicable
        ++mUmiTypeCounts[umiReadType.ordinal()];

        if(supportsVariant)
            ++mUmiTypeCounts[umiReadType.ordinal() + 3];
    }

    private void registerRawSupport(final RawContext rawContext, double recalibratedBaseQuality)
    {
        if(rawContext.AltSupport)
        {
            ++mRawAltSupport;
            mRawAltBaseQualityTotal += rawContext.BaseQuality;
            mRecalibratedBaseQualityTotal += recalibratedBaseQuality;
        }
        else if(rawContext.RefSupport)
        {
            ++mRawRefSupport;
            mRawRefBaseQualityTotal += rawContext.BaseQuality;
            mRecalibratedBaseQualityTotal += recalibratedBaseQuality;
        }
    }

    private void addVariantVisRecord(
            final SAMRecord record, final MatchType matchType,
            @Nullable QualityCalculator.QualityScores modifiedQualities, @Nullable final FragmentData fragmentData)
    {
        if(mVariantVis != null)
            mVariantVis.addEvidence(record, fragmentData, matchType, modifiedQualities);
    }

    private void logReadEvidence(final SAMRecord record, final MatchType matchType, int readIndex, double quality)
    {
        if(!mConfig.LogEvidenceReads || !SG_LOGGER.isTraceEnabled())
            return;

        // mQualityCalculator.logReadQualCalcs(this, readIndex, record, adjustedNumOfEvents);

        // Variant,MatchType,ReadId,ReadStart,Cigar,LeftCore,Index,RightCore,ReadIndex,Quality
        SG_LOGGER.trace("READ_EV,{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                mSample, chromosome(), position(), ref(), alt(),
                matchType, record.getReadName(), record.getAlignmentStart(), record.getCigarString(),
                mReadContext.indexedBases().LeftCoreIndex, mReadContext.indexedBases().Index,
                mReadContext.indexedBases().RightCoreIndex, readIndex, format("%.1f", quality));
    }

    private RawContext createRawContextFromCoreMatch(final SAMRecord record)
    {
        // check for an exact core match by which to centre the read
        if(mMaxCandidateDeleteLength < 5)
            return RawContext.INVALID_CONTEXT;

        int scLenLeft = leftSoftClipLength(record);
        int scLenRight = rightSoftClipLength(record);

        if(max(scLenLeft, scLenRight) < 5)
            return RawContext.INVALID_CONTEXT;

        final String variantCore = mReadContext.coreString();
        int coreStartIndex = record.getReadString().indexOf(variantCore);
        if(coreStartIndex < 1)
            return RawContext.INVALID_CONTEXT;

        int coreEndIndex = coreStartIndex + variantCore.length() - 1;
        boolean isValidRead = false;
        int baseQuality = 0;

        if(scLenLeft > scLenRight)
        {
            // the core match must span from the left soft-clipping into the matched bases
            int readRightCorePosition = record.getReferencePositionAtReadPosition(coreEndIndex + 1);

            if(readRightCorePosition > mVariant.position() && coreStartIndex < scLenLeft)
            {
                isValidRead = true;
                baseQuality = record.getBaseQualities()[scLenLeft];
            }
        }
        else
        {
            // the core match must span from the matched bases into the right soft-clipping
            int readLeftCorePosition = record.getReferencePositionAtReadPosition(coreStartIndex + 1);
            int postCoreIndexDiff = record.getReadBases().length - coreEndIndex;

            if(readLeftCorePosition > 0 && readLeftCorePosition < mVariant.position() && postCoreIndexDiff < scLenRight)
            {
                isValidRead = true;
                baseQuality = record.getBaseQualities()[record.getReadBases().length - scLenRight];
            }
        }

        if(!isValidRead)
            return RawContext.INVALID_CONTEXT;

        int readIndex = coreStartIndex + mReadContext.indexedBases().Index - mReadContext.indexedBases().LeftCoreIndex;

        return new RawContext(
                readIndex, false, false, true,
                true, false, true, baseQuality);
    }

    private boolean readHasIndelInCore(final SAMRecord record)
    {
        if(!record.getCigar().containsOperator(D) && !record.getCigar().containsOperator(I))
            return false;

        int variantLeftCorePos = mVariant.position() - (mReadContext.indexedBases().Index - mReadContext.indexedBases().LeftCoreIndex);
        int variantRightCorePos = mVariant.position() + (mReadContext.indexedBases().RightCoreIndex - mReadContext.indexedBases().Index);

        int currentPos = record.getAlignmentStart() - 1;

        // eg 2S10M2D10M starting at 100: first non-SC element, in this case a delete, starts at 109
        for(CigarElement element : record.getCigar())
        {
            if(element.getOperator() == S)
                continue;

            if(element.getOperator() == I || element.getOperator() == D)
            {
                int indelLowerPos = currentPos;
                int indelUpperPos = indelLowerPos + (element.getOperator() == D ? element.getLength() : 1);

                if(positionsOverlap(variantLeftCorePos, variantRightCorePos, indelLowerPos, indelUpperPos))
                    return true;
            }
            else if(element.getOperator() == M)
            {
                currentPos += element.getLength();
            }
        }

        return false;
    }

    private int calcLeftAlignmentIndex(final SAMRecord record)
    {
        // Left alignment: Match full read context starting at base = pos - rc_index
        int leftCoreOffset = mReadContext.indexedBases().Index - mReadContext.indexedBases().LeftCoreIndex;
        int realignLeftCorePos = position() - leftCoreOffset;
        int realignLeftCoreIndex = record.getReadPositionAtReferencePosition(realignLeftCorePos);

        if(realignLeftCoreIndex > 0)
        {
            int realignLeftReadIndex = realignLeftCoreIndex - 1 + leftCoreOffset;
            return realignLeftReadIndex;
        }

        int deleteCount = (int)record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == D).count();
        if(deleteCount == 0)
            return -1;

        int deleteStartPos = record.getAlignmentStart();
        int deleteStartIndex = 0;
        for(CigarElement element : record.getCigar())
        {
            if(element.getOperator() == S || element.getOperator() == I)
                continue;

            if(element.getOperator() == M)
            {
                deleteStartPos += element.getLength();
                deleteStartIndex += element.getLength();
            }
            else if(element.getOperator() == D)
            {
                --deleteCount;

                if(deleteCount == 0)
                    break;

                deleteStartPos += element.getLength();
            }
        }

        --deleteStartPos;

        int posDiff = realignLeftCorePos - deleteStartPos;
        int realignLeftReadIndex = deleteStartIndex + posDiff - 1 + leftCoreOffset;
        return realignLeftReadIndex;
    }

    private int calcRightAlignmentIndex(final SAMRecord record)
    {
        // Right alignment: Match full read context ending at base = pos + length[RC} - rc_index - 1 - length(alt) + length(ref)
        int rightCoreOffset = mReadContext.indexedBases().RightCoreIndex - mReadContext.indexedBases().Index;
        int realignRightPos = position() + rightCoreOffset - mVariant.alt().length() + mVariant.ref().length();
        int realignRightCoreIndex = record.getReadPositionAtReferencePosition(realignRightPos);

        if(realignRightCoreIndex > 0)
        {
            int realignRightReadIndex = realignRightCoreIndex - 1 - rightCoreOffset;
            return realignRightReadIndex;
        }

        return -1;
    }

    private RealignedContext checkRealignment(final SAMRecord record)
    {
        // try left and right alignment in turn
        int realignLeftReadIndex = calcLeftAlignmentIndex(record);

        if(realignLeftReadIndex >= 0) //  && realignLeftReadIndex != readIndex
        {
            ReadContextMatch match = determineReadContextMatch(record, realignLeftReadIndex, false);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL)
                return new RealignedContext(EXACT, mReadContext.indexedBases().length(), realignLeftReadIndex);
        }

        int realignRightReadIndex = calcRightAlignmentIndex(record);

        if(realignRightReadIndex >= 0)
        {
            // still need to test even if this index matches the original readIndex since if the readIndex was in a delete
            // it will be have skipped above
            ReadContextMatch match = determineReadContextMatch(record, realignRightReadIndex, false);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL)
                return new RealignedContext(RealignedType.EXACT, mReadContext.indexedBases().length(), realignRightReadIndex);
            else if(match == ReadContextMatch.CORE_PARTIAL)
                return new RealignedContext(RealignedType.CORE_PARTIAL, mReadContext.indexedBases().length(), realignRightReadIndex);
        }

        // try a simple string search and take it as exact if the matched index is within the expected range
        String readContext = mReadContext.indexedBases().fullString();
        if(readContext.length() >= REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH)
        {
            int matchedReadIndex = record.getReadString().indexOf(readContext);

            if(matchedReadIndex >= 0)
            {
                int matchedIndex = matchedReadIndex + mReadContext.indexedBases().Index - mReadContext.indexedBases().LeftFlankIndex;
                if(abs(matchedIndex - realignLeftReadIndex) <= REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER
                || abs(matchedIndex - realignRightReadIndex) <= REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER)
                {
                    return new RealignedContext(RealignedType.EXACT, mReadContext.indexedBases().length(), matchedIndex);
                }
            }
        }

        return RealignedContext.NONE;
    }

    public void addLocalPhaseSet(int lps, int readCount, double allocCount)
    {
        if(mLocalPhaseSets == null)
        {
            mLocalPhaseSets = Lists.newArrayList();
            mLpsCounts = Lists.newArrayList();
        }

        // add in order of highest counts
        int index = 0;
        while(index < mLpsCounts.size())
        {
            int existingCount = mLpsCounts.get(index);
            if(readCount + allocCount > existingCount)
                break;

            ++index;
        }

        mLocalPhaseSets.add(index, lps);
        int lpsTotalCount = readCount + (int)allocCount;
        mLpsCounts.add(index, lpsTotalCount);
    }

    private void countAltSupportMetrics(final SAMRecord record, final FragmentData fragmentData)
    {
        mAltReadStrandBias.registerRead(record, fragmentData, mVariant);
        mAltFragmentStrandBias.registerFragment(record);

        if(mFragmentCoords != null)
        {
            if(fragmentData != null)
                mFragmentCoords.addRead(fragmentData.First, fragmentData.Second);
            else
                mFragmentCoords.addRead(record, null);
        }
    }

    private int getMaxRealignDistance(final SAMRecord record)
    {
        int index = mReadContext.readBasesPositionIndex();
        int leftIndex = mReadContext.readBasesLeftCentreIndex();
        int rightIndex = mReadContext.readBasesRightCentreIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        int indelLength = record.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator() == I || x.getOperator() == D).mapToInt(x -> x.getLength()).sum();

        return Math.max(indelLength + Math.max(leftOffset, rightOffset), Realignment.MAX_REPEAT_SIZE) + 1;
    }

    private void checkImproperCount(final SAMRecord record)
    {
        if(isImproperPair(record) || record.getSupplementaryAlignmentFlag())
        {
            mImproperPairCount++;
        }
    }

    public void applyMapQualityRatio()
    {
        int depth = depth();
        int avgTotalMapQuality = depth > 0 ? (int)round(mapQualityTotal() / (double)depth) : 0;

        if(avgTotalMapQuality == 0)
            return;

        int altSupport = altSupport();
        int avgAltMapQuality = altSupport > 0 ? (int)round(altMapQualityTotal() / (double)altSupport) : 0;

        double ratioRaw = (avgAltMapQuality + MQ_RATIO_SMOOTHING) / (avgTotalMapQuality + MQ_RATIO_SMOOTHING);
        double calcRatio = pow(min(1, ratioRaw), mConfig.Quality.MapQualityRatioFactor);

        mQualities.applyRatio(calcRatio);
    }

    public boolean logEvidence() { return mConfig.LogEvidenceReads; }

    @VisibleForTesting
    public ReadSupportCounts readSupportQualityCounts() { return mQualities; };
    public ReadSupportCounts readSupportCounts() { return mCounts; }
    public FragmentCoords fragmentCoords() { return mFragmentCoords; }
}
