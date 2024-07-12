package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.CORE;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.FULL;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REALIGNED;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REF;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.EVIDENCE_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MQ_RATIO_SMOOTHING;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;
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
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;
import static com.hartwig.hmftools.sage.evidence.Realignment.INVALID_INDEX;
import static com.hartwig.hmftools.sage.evidence.Realignment.checkRealignment;
import static com.hartwig.hmftools.sage.evidence.Realignment.considerRealignedDel;
import static com.hartwig.hmftools.sage.evidence.Realignment.realignedReadIndexPosition;
import static com.hartwig.hmftools.sage.evidence.SplitReadSegment.formSegment;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.DELETED;
import static com.hartwig.hmftools.sage.filter.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.N;

import java.util.Arrays;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.filter.FragmentCoords;
import com.hartwig.hmftools.sage.filter.StrandBiasData;
import com.hartwig.hmftools.sage.quality.ArtefactContext;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.quality.ReadContextQualCache;
import com.hartwig.hmftools.sage.quality.UltimaQualModel;
import com.hartwig.hmftools.sage.common.NumberEvents;
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
    private final boolean mIsReferenceSample;
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

    private List<Integer> mLocalPhaseSets;
    private List<Integer> mLpsCounts;
    private int[] mUmiTypeCounts;
    private FragmentLengthCounts mFragmentLengthData;
    private FragmentCoords mFragmentCoords;

    // info only for VCF
    private double mTumorQualProbability;
    private double mMapQualFactor;

    public ReadContextCounter(
            final int id, final VariantReadContext readContext, final VariantTier tier, int maxCoverage, int minNumberOfEvents,
            final SageConfig config, final QualityCalculator qualityCalculator, final String sampleId, boolean isReferenceSample)
    {
        mId = id;

        mTier = tier;
        mMaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;
        mSample = sampleId;
        mIsReferenceSample = isReferenceSample;
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
        mUmiTypeCounts = null;
        mFragmentLengthData = mConfig.WriteFragmentLengths ? new FragmentLengthCounts() : null;
        mFragmentCoords = mConfig.Quality.HighDepthMode ? new FragmentCoords(REQUIRED_UNIQUE_FRAG_COORDS) : null;

        mTumorQualProbability = 0;
        mMapQualFactor = 0;
    }

    public int id() { return mId; }
    public SimpleVariant variant() { return mVariant; }
    public VariantReadContext readContext() { return mReadContext; }
    public ReadContextMatcher matcher() { return mMatcher; }
    public VariantTier tier() { return mTier; }
    public boolean isSnv() { return mVariant.isSNV(); }
    public boolean isIndel() { return mIsIndel; }
    public boolean isLongInsert() { return SimpleVariant.isLongInsert(mVariant); }
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

    public int depth() { return mCounts.Total; }

    public QualCounters qualCounters() { return mQualCounters; }

    public int baseQualityTotal() { return mQualCounters.baseQualityTotal(); }
    public int altBaseQualityTotal() { return mQualCounters.altBaseQualityTotal(); }

    public long mapQualityTotal() { return mQualCounters.mapQualityTotal(); }
    public long altMapQualityTotal() { return mQualCounters.altMapQualityTotal(); }

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

    public ArtefactContext artefactContext() { return mReadContext.artefactContext(); }
    public UltimaQualModel ultimaQualModel() { return mReadContext.ultimaQualModel(); }
    public boolean useMsiErrorRate() { return mQualCache.msiIndelErrorQual() != INVALID_BASE_QUAL;}

    public StrandBiasData fragmentStrandBiasAlt() { return mAltFragmentStrandBias; }
    public StrandBiasData fragmentStrandBiasNonAlt() { return mNonAltFragmentStrandBias; }
    public StrandBiasData readStrandBiasAlt() { return mAltReadStrandBias; }
    public StrandBiasData readStrandBiasNonAlt() { return mNonAltReadStrandBias; }
    public int improperPairCount() { return mImproperPairCount; }

    public ReadEdgeDistance readEdgeDistance() { return mReadEdgeDistance; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    @Nullable
    public VariantVis variantVis() { return mVariantVis; }

    public double averageAltBaseQuality()
    {
        // excludes realigned
        int supportCount = mCounts.Full + mCounts.PartialCore + mCounts.Core + mCounts.Realigned;
        return supportCount > 0 ? mQualCounters.altBaseQualityTotal() / (double)supportCount : 0;
    }

    public void setMaxCandidateDeleteLength(int length) { mMaxCandidateDeleteLength = length; }
    public int maxCandidateDeleteLength() { return mMaxCandidateDeleteLength; }
    public int maxPositionVsReadStart() { return mMaxPositionVsReadStart; }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<Integer> lpsCounts() { return mLpsCounts; }

    public int[] umiTypeCounts() { return mUmiTypeCounts; }
    public FragmentLengthCounts fragmentLengths() { return mFragmentLengthData; }

    public boolean exceedsMaxCoverage() { return mCounts.Total >= mMaxCoverage; }

    public boolean belowMinFragmentCoords() { return mFragmentCoords != null && !mFragmentCoords.atCapacity(); }

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

        boolean coreCovered = coversVariant(record, readVarIndex, splitReadSegment);

        if(!coreCovered && !checkRealigned)
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

        ReadContextMatch matchType = NONE;
        double calcBaseQuality = 0;
        double modifiedQuality = 0;
        QualityScores qualityScores = null;

        if(coreCovered)
        {
            calcBaseQuality = mQualityCalculator.calculateBaseQuality(this, readVarIndex, record);

            if(belowQualThreshold(calcBaseQuality))
            {
                if(rawContext.PositionType != DELETED)
                {
                    matchType = determineReadContextMatch(record, readVarIndex, splitReadSegment);

                    if(matchType.SupportsAlt)
                        countAltSupportMetrics(record, fragmentData);
                }

                addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                return UNRELATED;
            }

            qualityScores = mQualityCalculator.calculateQualityScores(
                    this, readVarIndex, record, adjustedNumOfEvents, calcBaseQuality);

            modifiedQuality = qualityScores.ModifiedQuality;

            // check for alt support (full, partial core or core)
            matchType = rawContext.PositionType != DELETED ? determineReadContextMatch(record, readVarIndex, splitReadSegment) : NONE;

            if(matchType.SupportsAlt)
            {
                VariantReadSupport readSupport = matchType.toReadSupport();

                registerReadSupport(record, readSupport, modifiedQuality, readVarIndex);

                mQualCounters.update(qualityScores, record.getMappingQuality(), matchType);

                mReadEdgeDistance.update(record, fragmentData, true);

                addVariantVisRecord(record, matchType, qualityScores, fragmentData);
                logReadEvidence(record, matchType, readVarIndex, modifiedQuality);
                countAltSupportMetrics(record, fragmentData);

                return ALT_SUPPORT;
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

                qualityScores = mQualityCalculator.calculateQualityScores(
                        this, realignedReadIndex, record, adjustedNumOfEvents, calcBaseQuality);

                modifiedQuality = qualityScores.ModifiedQuality;

                if(belowQualThreshold(calcBaseQuality))
                {
                    addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
                    return UNRELATED;
                }

                if(realignedType == EXACT)
                {
                    matchType = ReadContextMatch.REALIGNED;
                    registerReadSupport(record, REALIGNED, modifiedQuality, readVarIndex);

                    mQualCounters.update(qualityScores, record.getMappingQuality(), matchType);

                    addVariantVisRecord(record, matchType, qualityScores, fragmentData);
                    logReadEvidence(record, matchType, readVarIndex, modifiedQuality);

                    return ALT_SUPPORT;
                }

                if(realignedType == SHORTENED)
                    mJitterData.update(JitterMatch.SHORTENED);
                else if(realignedType == LENGTHENED)
                    mJitterData.update(JitterMatch.LENGTHENED);
            }
            else if(!coreCovered || readVarIndex < 0)
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

        registerReadSupport(record, readSupport, modifiedQuality, readVarIndex);
        mReadEdgeDistance.update(record, fragmentData, false);

        addVariantVisRecord(record, matchType, qualityScores, fragmentData);
        logReadEvidence(record, matchType, readVarIndex, modifiedQuality);

        return matchType == ReadContextMatch.REF ? REF_SUPPORT : UNRELATED;
    }

    private boolean coversVariant(final SAMRecord record, int readIndex, final SplitReadSegment splitReadSegment)
    {
        if(splitReadSegment != null)
            return mMatcher.coversVariant(splitReadSegment.ReadBases, splitReadSegment.ReadVarIndex);

        return mMatcher.coversVariant(record.getReadBases(), readIndex);
    }

    private ReadContextMatch determineReadContextMatch(final SAMRecord record, int readIndex, final SplitReadSegment splitReadSegment)
    {
        if(splitReadSegment != null)
        {
            return mMatcher.determineReadMatch(
                    splitReadSegment.ReadBases, splitReadSegment.ReadQuals, splitReadSegment.ReadVarIndex, false);
        }

        return mMatcher.determineReadMatch(record, readIndex);
    }

    private boolean belowQualThreshold(double calcBaseQuality)
    {
        return !mQualCache.usesMsiIndelErrorQual() && mConfig.Quality.HighDepthMode && calcBaseQuality < mConfig.Quality.HighBaseQualLimit;
    }

    private void registerReadSupport(final SAMRecord record, @Nullable final VariantReadSupport support, double quality, int readVarIndex)
    {
        mCounts.addSupport(support, 1);
        mQualities.addSupport(support, (int) quality);

        boolean supportsVariant = support != null
                && (support == FULL || support == VariantReadSupport.PARTIAL_CORE || support == CORE || support == REALIGNED);

        if(mConfig.Sequencing.HasUMIs)
        {
            countUmiType(record, supportsVariant);
        }

        if(mFragmentLengthData != null && (support == REF || supportsVariant))
        {
            mFragmentLengthData.addLength(abs(record.getInferredInsertSize()), supportsVariant);
        }

        if(support != null && (support == FULL || support == VariantReadSupport.PARTIAL_CORE))
        {
            if(isImproperPair(record) || record.getSupplementaryAlignmentFlag())
                mImproperPairCount++;
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

    public double mapQualFactor() { return mMapQualFactor; }
    public void setMapQualFactor(double factor) { mMapQualFactor = factor; }

    @VisibleForTesting
    public ReadSupportCounts readSupportQualityCounts() { return mQualities; };
    public ReadSupportCounts readSupportCounts() { return mCounts; }
    public FragmentCoords fragmentCoords() { return mFragmentCoords; }
}
