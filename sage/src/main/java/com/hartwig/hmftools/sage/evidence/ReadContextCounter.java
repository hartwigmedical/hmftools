package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.variant.SageVcfTags.CONSENSUS_TAG_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.CONSENSUS_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.CORE;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.FULL;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REALIGNED;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REF;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.isSbx;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.EVIDENCE_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.LONG_INSERT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.LONG_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MQ_RATIO_SMOOTHING;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS_2;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.SageConstants.TQP_QUAL_LOG_MIN;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatcher.isSimpleAltMatch;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;
import static com.hartwig.hmftools.sage.evidence.ReadEdgeDistance.calcAdjustedVariantPosition;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.ALT_SUPPORT_EXACT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.ALT_SUPPORT_LOW_QUAL_MISMATCHES;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.CHIMERIC;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.IN_SPLIT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAP_QUAL;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAX_COVERAGE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NON_CORE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_ALT_REF_MATCH;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.REF_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SOFT_CLIP;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.UNRELATED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LOW_QUAL_MISMATCHES;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;
import static com.hartwig.hmftools.sage.evidence.Realignment.INVALID_INDEX;
import static com.hartwig.hmftools.sage.evidence.Realignment.checkRealignment;
import static com.hartwig.hmftools.sage.evidence.Realignment.considerRealignedDel;
import static com.hartwig.hmftools.sage.evidence.Realignment.realignedReadIndexPosition;
import static com.hartwig.hmftools.sage.evidence.SplitReadSegment.formSegment;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.DELETED;
import static com.hartwig.hmftools.sage.filter.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isHighBaseQual;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.N;

import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.ReadMatchInfo;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.filter.FragmentCoords;
import com.hartwig.hmftools.sage.filter.FragmentLengths;
import com.hartwig.hmftools.sage.filter.StrandBiasData;
import com.hartwig.hmftools.sage.seqtech.IlluminaArtefactContext;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.quality.ReadContextQualCache;
import com.hartwig.hmftools.sage.common.NumberEvents;
import com.hartwig.hmftools.sage.seqtech.UltimaRealignedQualModels;
import com.hartwig.hmftools.sage.seqtech.UltimaVariantData;
import com.hartwig.hmftools.sage.vis.VariantVis;
import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounter
{
    private final int mId;
    private final VariantTier mTier;
    private final SimpleVariant mVariant;
    private final VariantReadContext mReadContext;
    private final ReadContextMatcher mMatcher;
    private final SageConfig mConfig;
    private final QualityCalculator mQualityCalculator;
    private final String mSample;
    private final int mMaxCoverage;
    private final VariantVis mVariantVis;

    // local variant-related state
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;
    private final boolean mIsIndel;
    private final ReadContextQualCache mQualCache;

    // counts and quals by support type
    private final ReadSupportCounts mQualities;
    private final ReadSupportCounts mCounts;
    private int mSimpleAltMatches;
    private int mHighQualStrongSupport;

    private final StrandBiasData mAltFragmentStrandBias;
    private final StrandBiasData mNonAltFragmentStrandBias;
    private final StrandBiasData mAltReadStrandBias;
    private final StrandBiasData mNonAltReadStrandBias;

    private final JitterData mJitterData;
    private int mImproperPairCount;

    private final QualCounters mQualCounters;

    private int mMaxCandidateDeleteLength;
    private final ReadEdgeDistance mReadEdgeDistance;
    private final int mMaxPositionVsReadStart;
    private int mNonAltNmCountTotal;
    private List<Integer> mLocalPhaseSets;
    private List<Integer> mLpsCounts;
    private int[] mConsensusTypeCounts;
    private FragmentLengthCounts mFragmentLengthData;
    private FragmentCoords mFragmentCoords;
    private final FragmentLengths mFragmentLengths;
    private double mAdjustedRefVaf;

    // info only for VCF
    private double mTumorQualProbability;
    private double mMapQualFactor;

    private final UltimaVariantData mUltimaData;

    public ReadContextCounter(
            final int id, final VariantReadContext readContext, final VariantTier tier, int maxCoverage, int minNumberOfEvents,
            final SageConfig config, final QualityCalculator qualityCalculator, final String sampleId, boolean isReferenceSample)
    {
        mId = id;

        mTier = tier;
        mMaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;
        mSample = sampleId;
        mQualityCalculator = qualityCalculator;
        mConfig = config;

        mReadContext = readContext;
        mMatcher = new ReadContextMatcher(mReadContext, true, isReferenceSample);
        mVariant = readContext.variant();

        mVariantVis = config.Visualiser.Enabled && config.Visualiser.processVariant(mVariant)
                ? new VariantVis(mConfig, mSample, mVariant, mReadContext, mTier) : null;

        // set local state to avoid testing on each read
        mIsMnv = mVariant.isMNV();
        mIsIndel = mVariant.isIndel();
        mQualCache = new ReadContextQualCache(readContext, qualityCalculator, sampleId);

        mQualities = new ReadSupportCounts();
        mCounts = new ReadSupportCounts();
        mSimpleAltMatches = 0;
        mHighQualStrongSupport = 0;

        mJitterData = new JitterData();

        mAltFragmentStrandBias = new StrandBiasData(true);
        mNonAltFragmentStrandBias = new StrandBiasData(false);
        mAltReadStrandBias = new StrandBiasData(true);
        mNonAltReadStrandBias = new StrandBiasData(false);

        mImproperPairCount = 0;

        mQualCounters = new QualCounters();
        mMaxCandidateDeleteLength = 0;
        mMaxPositionVsReadStart = mVariant.isDelete() ? mVariant.Position + abs(mVariant.indelLength()) : mVariant.Position;

        mReadEdgeDistance = new ReadEdgeDistance(calcAdjustedVariantPosition(mVariant.position(), variant().indelLengthAbs() + 1));

        mLocalPhaseSets = null;
        mLpsCounts = null;
        mConsensusTypeCounts = null;
        mFragmentLengthData = mConfig.WriteFragmentLengths ? new FragmentLengthCounts() : null;
        mFragmentCoords = new FragmentCoords(REQUIRED_UNIQUE_FRAG_COORDS_2);
        mFragmentLengths = new FragmentLengths();
        mAdjustedRefVaf = 0;

        mTumorQualProbability = 0;
        mMapQualFactor = 0;

        mUltimaData = isUltima() ? new UltimaVariantData(mReadContext) : null;
    }

    public int id() { return mId; }
    public SimpleVariant variant() { return mVariant; }
    public VariantReadContext readContext() { return mReadContext; }
    public ReadContextMatcher matcher() { return mMatcher; }
    public VariantTier tier() { return mTier; }
    public boolean isSnv() { return mVariant.isSNV(); }
    public boolean isIndel() { return mIsIndel; }
    public boolean isLongInsert() { return SageVariant.isLongInsert(mVariant); }
    public boolean isLongIndel() { return mVariant.indelLengthAbs() >= LONG_INSERT_LENGTH; }
    public boolean isInLongRepeat() { return mReadContext.maxRepeatCount() >= LONG_REPEAT_LENGTH; }
    public final ReadContextQualCache qualCache() { return mQualCache; }
    public String chromosome() { return mVariant.chromosome(); }
    public int position() { return mVariant.position(); }
    public String ref() { return mVariant.ref(); }
    public String alt() { return mVariant.alt(); }
    public String sampleId() { return mSample; }

    public int altSupport() { return mCounts.altSupport(); }
    public int strongAltSupport() { return mCounts.strongSupport(); }
    public int refSupport() { return mCounts.Ref; }

    public int simpleAltMatches() { return mSimpleAltMatches; }
    public int strongHighQualSupport() { return mHighQualStrongSupport; }

    public int depth() { return mCounts.Total; }

    public QualCounters qualCounters() { return mQualCounters; }

    public int baseQualityTotal() { return mQualCounters.baseQualityTotal(); }
    public int altBaseQualityTotal() { return mQualCounters.altRecalibratedBaseQualityTotal(); }

    public long mapQualityTotal() { return mQualCounters.mapQualityTotal(); }
    public long altMapQualityTotal() { return mQualCounters.altMapQualityTotal(); }
    public long nonAltNmCountTotal() { return mNonAltNmCountTotal; }
    public double vaf()
    {
        return mCounts.Total == 0 ? 0d : mCounts.altSupport() / (double)mCounts.Total;
    }

    public double tumorQuality()
    {
        int qualTotal = mQualities.Full + mQualities.PartialCore + mQualities.Realigned;

        if(!mJitterData.filterOnNoise())
            qualTotal *= mJitterData.qualBoost();

        return qualTotal;
    }

    public int[] counts() { return mCounts.toArray(); }
    public int[] quality() { return mQualities.toArray(); }
    public ReadSupportCounts readCounts() { return mCounts; }
    public ReadSupportCounts readQuals() { return mQualities; }

    public JitterData jitter() { return mJitterData; }

    public IlluminaArtefactContext artefactContext() { return mReadContext.artefactContext(); }
    public boolean useMsiErrorRate() { return mQualCache.msiIndelErrorQual() != INVALID_BASE_QUAL;}

    public StrandBiasData fragmentStrandBiasAlt() { return mAltFragmentStrandBias; }
    public StrandBiasData fragmentStrandBiasNonAlt() { return mNonAltFragmentStrandBias; }
    public StrandBiasData readStrandBiasAlt() { return mAltReadStrandBias; }
    public StrandBiasData readStrandBiasNonAlt() { return mNonAltReadStrandBias; }
    public FragmentLengths fragmentLengths() { return mFragmentLengths; }

    public int improperPairCount() { return mImproperPairCount; }

    public ReadEdgeDistance readEdgeDistance() { return mReadEdgeDistance; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    @Nullable
    public VariantVis variantVis() { return mVariantVis; }

    public double averageAltBaseQuality()
    {
        // excludes realigned
        int supportCount = mCounts.Full + mCounts.PartialCore + mCounts.Core + mCounts.Realigned + mQualCounters.lowQualAltSupportCount();
        return supportCount > 0 ? mQualCounters.altBaseQualityTotal() / (double)supportCount : 0;
    }

    public double averageAltRecalibratedBaseQuality()
    {
        // excludes realigned
        int supportCount = mCounts.Full + mCounts.PartialCore + mCounts.Core + mCounts.Realigned;
        return supportCount > 0 ? mQualCounters.altRecalibratedBaseQualityTotal() / (double)supportCount : 0;
    }

    public void setMaxCandidateDeleteLength(int length) { mMaxCandidateDeleteLength = length; }
    public int maxCandidateDeleteLength() { return mMaxCandidateDeleteLength; }
    public int maxPositionVsReadStart() { return mMaxPositionVsReadStart; }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<Integer> lpsCounts() { return mLpsCounts; }

    public int[] consensusTypeCounts() { return mConsensusTypeCounts; }
    public FragmentLengthCounts fragmentLengthCounts() { return mFragmentLengthData; }

    public UltimaVariantData ultimaData() { return mUltimaData; }

    public boolean exceedsMaxCoverage() { return mCounts.Total >= mMaxCoverage; }

    public String toString()
    {
        return format("id(%d) var(%s) core(%s) counts(f=%d p=%d c=%d)",
                mId, varString(), mReadContext.toString(), mCounts.Full, mCounts.PartialCore, mCounts.Core);
    }

    public String varString()
    {
        return format("%s:%d %s>%s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
    }

    public ReadMatchType processRead(final SAMRecord record, int numberOfEvents, @Nullable final FragmentData fragmentData)
    {
        if(exceedsMaxCoverage())
            return MAX_COVERAGE;

        if(mTier != VariantTier.HOTSPOT && record.getMappingQuality() < EVIDENCE_MIN_MAP_QUAL)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return MAP_QUAL;
        }

        if(mConfig.Quality.HighDepthMode && fragmentData == null && isChimericRead(record))
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return CHIMERIC;
        }

        RawContext rawContext = RawContext.createFromRead(mVariant, record);

        int readVarIndex = rawContext.ReadVariantIndex;
        SplitReadSegment splitReadSegment = record.getCigar().containsOperator(N) ? formSegment(record, mVariant.Position, readVarIndex) : null;

        if(rawContext.PositionType == VariantReadPositionType.LOW_QUAL)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return SOFT_CLIP;
        }

        boolean checkRealigned = false;
        Integer realignedReadIndex = null;

        if(rawContext.ReadVariantIndex < 0)
        {
            if(mVariant.isDelete() && considerRealignedDel(record, mMaxPositionVsReadStart))
            {
                realignedReadIndex = realignedReadIndexPosition(mReadContext, record);
                checkRealigned = realignedReadIndex != INVALID_INDEX;
            }

            if(!checkRealigned)
            {
                addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                return UNRELATED;
            }
        }

        if(rawContext.PositionType == VariantReadPositionType.SKIPPED)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return IN_SPLIT;
        }

        boolean variantCovered = coversVariant(record, readVarIndex, splitReadSegment);

        if(!variantCovered && !checkRealigned)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
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

        ReadMatchInfo readMatchInfo = ReadMatchInfo.NO_MATCH;
        ReadContextMatch matchType = NONE;
        double calcBaseQuality = 0;
        double modifiedQuality = 0;
        QualityScores qualityScores = null;

        if(variantCovered)
        {
            calcBaseQuality = mQualityCalculator.calculateBaseQuality(this, readVarIndex, record);

            // CHECK: what is the meaning of neg / invalid qual, what scenarios, why only added in Ultima - see check below too
            if(calcBaseQuality < 0)
            {
                addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                return UNRELATED;
            }

            qualityScores = mQualityCalculator.calculateQualityScores(
                    this, readVarIndex, record, adjustedNumOfEvents, calcBaseQuality);

            if(belowQualThreshold(calcBaseQuality))
            {
                if(rawContext.PositionType != DELETED)
                {
                    readMatchInfo = determineReadContextMatch(record, readVarIndex, splitReadSegment);
                    matchType = readMatchInfo.MatchType;

                    if(matchType.SupportsAlt)
                    {
                        countAltSupportMetrics(record, fragmentData);
                        mQualCounters.update(qualityScores);
                    }
                }

                addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                return UNRELATED;
            }

            modifiedQuality = qualityScores.ModifiedQuality;

            // check for alt support (full, partial core or core)
            if(rawContext.PositionType != DELETED)
            {
                readMatchInfo = determineReadContextMatch(record, readVarIndex, splitReadSegment);
                matchType = readMatchInfo.MatchType;
            }

            if(matchType.SupportsAlt)
            {
                VariantReadSupport readSupport = matchType.toReadSupport();

                registerReadSupport(record, readSupport, readVarIndex, modifiedQuality, calcBaseQuality);

                mQualCounters.update(qualityScores, record.getMappingQuality(), matchType);

                mReadEdgeDistance.update(record, fragmentData, matchType.FullAltSupport);

                addVariantVisRecord(record, matchType, qualityScores, fragmentData);
                logReadEvidence(record, matchType, readVarIndex, modifiedQuality);
                countAltSupportMetrics(record, fragmentData);

                return readMatchInfo.ExactMatch && matchType.FullAltSupport ? ALT_SUPPORT_EXACT : ALT_SUPPORT_LOW_QUAL_MISMATCHES;
            }
        }

        RealignedType realignedType = RealignedType.NONE;

        if(matchType != ReadContextMatch.REF)
        {
            if(realignedReadIndex == null)
                realignedReadIndex = realignedReadIndexPosition(mReadContext, record);

            boolean canRealign = realignedReadIndex != INVALID_INDEX && coversVariant(record, realignedReadIndex, splitReadSegment);

            if(canRealign)
                realignedType = checkRealignment(mReadContext, mMatcher, record, readVarIndex, realignedReadIndex, splitReadSegment);

            if(realignedType != RealignedType.NONE)
            {
                // recompute qual off this realigned index
                calcBaseQuality = mQualityCalculator.calculateBaseQuality(this, realignedReadIndex, record);

                if(calcBaseQuality < 0)
                {
                    addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                    return UNRELATED;
                }

                qualityScores = mQualityCalculator.calculateQualityScores(
                        this, realignedReadIndex, record, adjustedNumOfEvents, calcBaseQuality);

                modifiedQuality = qualityScores.ModifiedQuality;

                if(belowQualThreshold(calcBaseQuality))
                {
                    addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                    return UNRELATED;
                }

                if(realignedType == EXACT || realignedType == LOW_QUAL_MISMATCHES)
                {
                    matchType = ReadContextMatch.REALIGNED;
                    registerReadSupport(record, REALIGNED, readVarIndex, modifiedQuality, calcBaseQuality);

                    mQualCounters.update(qualityScores, record.getMappingQuality(), matchType);

                    mReadEdgeDistance.update(record, fragmentData, true);

                    addVariantVisRecord(record, matchType, qualityScores, fragmentData);
                    logReadEvidence(record, matchType, readVarIndex, modifiedQuality);
                    countAltSupportMetrics(record, fragmentData);

                    return realignedType == EXACT ? ALT_SUPPORT_EXACT : ALT_SUPPORT_LOW_QUAL_MISMATCHES;
                }

                if(realignedType == SHORTENED)
                    mJitterData.update(JitterMatch.SHORTENED);
                else if(realignedType == LENGTHENED)
                    mJitterData.update(JitterMatch.LENGTHENED);
            }
            else if(!variantCovered || readVarIndex < 0)
            {
                // exit if would have earlier but for the realignment test
                addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                return readVarIndex < 0 ? UNRELATED : NON_CORE;
            }
        }

        if(realignedType == RealignedType.NONE)
        {
            JitterMatch jitterMatch = checkJitter(mReadContext, mMatcher, record, readVarIndex);
            mJitterData.update(jitterMatch);
        }

        VariantReadSupport readSupport = null;

        if(matchType == ReadContextMatch.REF)
        {
            readSupport = REF;
        }
        else if(matchType == ReadContextMatch.SIMPLE_ALT && rawContext.PositionType != VariantReadPositionType.SOFT_CLIP)
        {
            ++mSimpleAltMatches;
        }

        // special case to ignore updating depth and other metrics when an indel could not have support the alt
        if(matchType == NONE && mVariant.isInsert() && readVarIndex - mVariant.indelLength() < 0)
            return UNRELATED;

        mQualCounters.update(qualityScores, record.getMappingQuality(), matchType);

        mNonAltFragmentStrandBias.registerFragment(record);
        mNonAltReadStrandBias.registerRead(record, fragmentData, this);

        registerReadSupport(record, readSupport, readVarIndex, modifiedQuality, calcBaseQuality);
        mReadEdgeDistance.update(record, fragmentData, false);
        mFragmentLengths.processRead(record, false);

        addVariantVisRecord(record, matchType, qualityScores, fragmentData);
        logReadEvidence(record, matchType, readVarIndex, modifiedQuality);
        mNonAltNmCountTotal += numberOfEvents;

        if(matchType == ReadContextMatch.REF)
            return REF_SUPPORT;

        // test for a simple non-match with the alt, for negative phasing
        Boolean simpleMatch = isSimpleAltMatch(mVariant, record, readVarIndex);

        if(simpleMatch == Boolean.FALSE)
            return NO_ALT_REF_MATCH;

        return UNRELATED;
    }

    private boolean coversVariant(final SAMRecord record, int readIndex, final SplitReadSegment splitReadSegment)
    {
        if(splitReadSegment != null)
            return mMatcher.coversVariant(splitReadSegment.ReadBases, splitReadSegment.ReadVarIndex);

        return mMatcher.coversVariant(record.getReadBases(), readIndex);
    }

    private ReadMatchInfo determineReadContextMatch(final SAMRecord record, int readIndex, final SplitReadSegment splitReadSegment)
    {
        if(splitReadSegment != null)
        {
            return mMatcher.determineReadMatchInfo(
                    splitReadSegment.ReadBases, splitReadSegment.ReadQuals, splitReadSegment.ReadVarIndex, false);
        }

        return mMatcher.determineReadMatchInfo(record, readIndex);
    }

    private boolean belowQualThreshold(double calcBaseQuality)
    {
        return mConfig.Quality.HighDepthMode && !mQualCache.usesMsiIndelErrorQual() && calcBaseQuality < mConfig.Quality.HighBaseQualLimit;
    }

    private void registerReadSupport(
            final SAMRecord record, @Nullable final VariantReadSupport support, int readVarIndex, double quality, double baseQuality)
    {
        boolean supportsAlt = false;
        boolean supportsAltStrong = false;

        if(support != null)
        {
            supportsAltStrong = support == FULL || support == VariantReadSupport.PARTIAL_CORE || support == REALIGNED;
            supportsAlt = supportsAltStrong || support == CORE;
        }

        // update UMI counts prior to read counts so if they need to be initialised, they do not double-count the current read
        countConsensusType(record, readVarIndex, supportsAltStrong);

        mCounts.addSupport(support, 1);
        mQualities.addSupport(support, (int)quality);

        if(supportsAltStrong && isHighBaseQual(baseQuality))
            ++mHighQualStrongSupport;

        if(mFragmentLengthData != null && (support == REF || supportsAlt))
        {
            mFragmentLengthData.addLength(abs(record.getInferredInsertSize()), supportsAlt);
        }

        if(support != null && (support == FULL || support == VariantReadSupport.PARTIAL_CORE))
        {
            if(isImproperPair(record) || record.getSupplementaryAlignmentFlag())
                mImproperPairCount++;
        }

        if(mUltimaData != null && support != null && support == FULL)
        {
            mUltimaData.addReadSupportInfo(record, readVarIndex);
        }
    }

    private void countConsensusType(final SAMRecord record, int readVarIndex, boolean supportsVariant)
    {
        if(mConsensusTypeCounts == null)
        {
            if(!record.hasAttribute(CONSENSUS_TYPE_ATTRIBUTE))
                return;

            // 3 total depth values followed by the 3 variant support values
            mConsensusTypeCounts = new int[CONSENSUS_TAG_TYPE_COUNT];

            mConsensusTypeCounts[ConsensusType.NONE.ordinal()] = depth();
            mConsensusTypeCounts[ConsensusType.NONE.ordinal() + CONSENSUS_TYPE_COUNT] = strongAltSupport();
        }

        ConsensusType consensusType = extractConsensusType(record);

        if(consensusType == DUAL && isSbx())
        {
            int duplexBaseIndex = SbxBamUtils.extractDuplexBaseIndex(record);

            if(!SbxBamUtils.inDuplexRegion(!record.getReadNegativeStrandFlag(), duplexBaseIndex, readVarIndex))
                consensusType = SINGLE;
        }

        // add to total and variant support if applicable
        ++mConsensusTypeCounts[consensusType.ordinal()];

        if(supportsVariant)
            ++mConsensusTypeCounts[consensusType.ordinal() + 3];
    }

    private void addVariantVisRecord(
            final SAMRecord record, final ReadContextMatch matchType,
            @Nullable QualityScores modifiedQualities, @Nullable final FragmentData fragmentData)
    {
        if(mVariantVis != null)
            mVariantVis.addEvidence(record, fragmentData, matchType, modifiedQualities);
    }

    private void logReadEvidence(final SAMRecord record, final ReadContextMatch matchType, int readIndex, double quality)
    {
        if(!mConfig.LogEvidenceReads || !SG_LOGGER.isTraceEnabled())
            return;

        // mQualityCalculator.logReadQualCalcs(this, readIndex, record, adjustedNumOfEvents);

        // Variant,MatchType,ReadId,ReadStart,Cigar,LeftCore,Index,RightCore,ReadIndex,Quality
        SG_LOGGER.trace("READ_EV,{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                mSample, chromosome(), position(), ref(), alt(),
                matchType, record.getReadName(), record.getAlignmentStart(), record.getCigarString(),
                mReadContext.CoreIndexStart, mReadContext.VarIndex,
                mReadContext.CoreIndexEnd, readIndex, format("%.1f", quality));
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
        mAltReadStrandBias.registerRead(record, fragmentData, this);
        mAltFragmentStrandBias.registerFragment(record);

        if(mFragmentCoords != null)
        {
            if(fragmentData != null)
                mFragmentCoords.addRead(fragmentData.First, fragmentData.Second);
            else
                mFragmentCoords.addRead(record, null);
        }

        mFragmentLengths.processRead(record, true);
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

    public double tumorQualProbability() { return mTumorQualProbability; }
    public void setTumorQualProbability(double probability) { mTumorQualProbability = probability; }
    public double logTqp()
    {
        double tqp = min(max(tumorQualProbability(), 0), 1);
        return log10(max(tqp, TQP_QUAL_LOG_MIN));
    }

    public double mapQualFactor() { return mMapQualFactor; }
    public void setMapQualFactor(double factor) { mMapQualFactor = factor; }

    public double adjustedRefVaf() { return mAdjustedRefVaf; }
    public void setAdjustedRefVaf(double vaf) { mAdjustedRefVaf = vaf; }
    public QualityCalculator qualityCalculator() { return mQualityCalculator; }

    @VisibleForTesting
    public ReadSupportCounts readSupportQualityCounts() { return mQualities; };
    public ReadSupportCounts readSupportCounts() { return mCounts; }
    public FragmentCoords fragmentCoords() { return mFragmentCoords; }
}
