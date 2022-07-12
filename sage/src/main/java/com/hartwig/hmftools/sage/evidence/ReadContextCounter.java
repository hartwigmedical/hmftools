package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_BASE_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_EVIDENCE_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.UNRELATED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.jitterPenalty;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.quality.QualityConfig;
import com.hartwig.hmftools.sage.SageConfig;
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
    public final VariantTier Tier;
    public final int MaxCoverage;

    private final int mId;

    private final VariantHotspot mVariant;
    private final ReadContext mReadContext;
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;

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

    private int mSoftClipInsertSupport;
    private int mMaxCandidateDeleteLength;

    private List<Integer> mLocalPhaseSets;
    private List<int[]> mLpsCounts;

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
            final int maxCoverage, final int minNumberOfEvents)
    {
        mId = id;

        Tier = tier;
        MaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;

        mReadContext = readContext;
        mVariant = variant;
        mIsMnv = variant.isMNV();

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
        mSoftClipInsertSupport = 0;
        mMaxCandidateDeleteLength = 0;

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

    public double vaf() { return alleleFrequency(altSupport()); }

    private double alleleFrequency(double support)
    {
        return mCounts[RC_TOTAL] == 0 ? 0d : support / mCounts[RC_TOTAL];
    }

    public int tumorQuality()
    {
        int tumorQuality = mQualities[RC_FULL] + mQualities[RC_PARTIAL] + mQualities[RC_REALIGNED];
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

    public int strandDepth() { return mForwardStrand + mReverseStrand; }

    public int improperPair() { return mImproperPair; }

    public int rawDepth() { return mRawDepth; }
    public int rawAltSupport() { return mRawAltSupport; }
    public int rawRefSupport() { return mRawRefSupport; }
    public int rawAltBaseQuality() { return mRawAltBaseQuality; }
    public int rawRefBaseQuality() { return mRawRefBaseQuality; }
    public int softClipInsertSupport() { return mSoftClipInsertSupport; }
    public void setMaxCandidateDeleteLength(int length) { mMaxCandidateDeleteLength = length; }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<int[]> lpsCounts() { return mLpsCounts; }

    public boolean exceedsMaxCoverage() { return mCounts[RC_TOTAL] >= MaxCoverage; }

    public String toString()
    {
        return format("id(%d) var(%s) core(%s) counts(f=%d p=%d c=%d)",
                mId, varString(), mReadContext.toString(), mCounts[RC_FULL], mCounts[RC_PARTIAL], mCounts[RC_CORE]);
    }

    private String varString()
    {
        return format("%s:%d %s>%s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
    }

    public ReadMatchType processRead(
            final SAMRecord record, final SageConfig sageConfig, final QualityCalculator qualityCalc, int numberOfEvents)
    {
        if(exceedsMaxCoverage())
            return UNRELATED;

        if(!Tier.equals(VariantTier.HOTSPOT) && record.getMappingQuality() < DEFAULT_EVIDENCE_MAP_QUAL)
            return UNRELATED;

        RawContext rawContext = RawContext.create(mVariant, record);

        if(rawContext.ReadIndex < 0)
        {
            if(rawContext.DepthSupport || rawContext.AltSupport || rawContext.RefSupport)
            {
                SG_LOGGER.error("rawContext missing readIndex but with support(depth={} ref={} alt={}",
                        rawContext.DepthSupport, rawContext.AltSupport, rawContext.RefSupport);
            }

            rawContext = createRawContextFromCoreMatch(record);

            if(rawContext.ReadIndex < 0)
                return UNRELATED;
        }

        if(rawContext.ReadIndexInSkipped)
            return UNRELATED;

        int readIndex = rawContext.ReadIndex;
        boolean baseDeleted = rawContext.ReadIndexInDelete;

        mRawDepth += rawContext.DepthSupport ? 1 : 0;
        mRawAltSupport += rawContext.AltSupport ? 1 : 0;
        mRawRefSupport += rawContext.RefSupport ? 1 : 0;
        mRawAltBaseQuality += rawContext.AltQuality;
        mRawRefBaseQuality += rawContext.RefQuality;

        if(rawContext.ReadIndexInSoftClip && rawContext.AltSupport)
            ++mSoftClipInsertSupport;

        boolean covered = mReadContext.indexedBases().isCoreCovered(readIndex, record.getReadBases().length);

        if(!covered)
            return UNRELATED;

        final QualityConfig qualityConfig = sageConfig.Quality;

        double adjustedNumOfEvents = numberOfEvents;

        if(mIsMnv)
            adjustedNumOfEvents = NumberEvents.calcWithMnvRaw(numberOfEvents, mVariant.ref(), mVariant.alt());

        if(max(mVariant.ref().length(), mVariant.alt().length()) <= SC_READ_EVENTS_FACTOR)
        {
            // penalise variants except long INDELs for their soft-clipped bases
            adjustedNumOfEvents += NumberEvents.calcSoftClipAdjustment(record);
        }

        adjustedNumOfEvents = max(mMinNumberOfEvents, adjustedNumOfEvents);

        double quality = qualityCalc.calculateQualityScore(this, readIndex, record, adjustedNumOfEvents);

        // Check if FULL, PARTIAL, OR CORE
        if(!baseDeleted)
        {
            boolean wildcardMatchInCore = mVariant.isSNV() && mReadContext.microhomology().isEmpty();

            int maxCoreMismatches = mVariant.isIndel() && mVariant.alt().length() >= CORE_LOW_QUAL_MISMATCH_BASE_LENGTH ?
                    mVariant.alt().length() / CORE_LOW_QUAL_MISMATCH_BASE_LENGTH : 0;

            IndexedBases readBases = record.getCigar().containsOperator(CigarOperator.N) ?
                    ExpandedBasesFactory.expand(position(), readIndex, record) :
                    new IndexedBases(position(), readIndex, record.getReadBases());

            final ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                    readBases, record.getBaseQualities(), wildcardMatchInCore, maxCoreMismatches);

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

                /*
                SG_LOGGER.trace("var({}) readContext({}-{}-{}) support({}) read(idx={} posStart={} cigar={} id={}) readBases({})",
                        varString(), mReadContext.indexedBases().LeftCoreIndex, mReadContext.indexedBases().Index,
                        mReadContext.indexedBases().RightCoreIndex, match, readIndex, record.getAlignmentStart(), record.getCigarString(),
                        record.getReadName(), record.getReadString());
                */

                countStrandedness(record);
                checkImproperCount(record);
                return SUPPORT;
            }
        }

        boolean canRealign = abs(mVariant.indelLength()) >= REALIGN_READ_MIN_INDEL_LENGTH || readHasIndelInCore(record);
        RealignedContext realignment = canRealign ? checkRealignment(record) : RealignedContext.NONE;

        /*
        // log differences
        RealignedContext oldRealignment = Realignment.realignedAroundIndex(mReadContext, readIndex, record.getReadBases(), getMaxRealignDistance(record));
        // if((realignment.Type == EXACT) != (oldRealignment.Type == EXACT))
        if(realignment.Type == NONE && oldRealignment.Type == EXACT && rawContext.ReadIndexInDelete)
        {
            SG_LOGGER.info("var({}) realign diff: new({}) old({} m={} mi={}) readContext({}-{}-{}) read(idx={} posStart={} cigar={} id={}) readCxt({}) readBases({})",
                    varString(), realignment.Type, oldRealignment.Type, oldRealignment.MatchLength, oldRealignment.MatchReadIndex,
                    mReadContext.indexedBases().LeftCoreIndex, mReadContext.indexedBases().Index,
                    mReadContext.indexedBases().RightCoreIndex,  readIndex, record.getAlignmentStart(), record.getCigarString(),
                    record.getReadName(), mReadContext.indexedBases().fullString(), record.getReadString());
        }
        */

        if(realignment.Type == EXACT)
        {
            mCounts[RC_REALIGNED]++;
            mQualities[RC_REALIGNED] += quality;

            mCounts[RC_TOTAL]++;
            mQualities[RC_TOTAL] += quality;
            return SUPPORT;
        }

        // switch back to the old method to test for jitter
        RealignedContext jitterRealign = canRealign ? Realignment.realignedAroundIndex(
                mReadContext, readIndex, record.getReadBases(), getMaxRealignDistance(record)) : RealignedContext.NONE;

        if(rawContext.ReadIndexInSoftClip && !rawContext.AltSupport)
        {
            if(jitterRealign.Type != LENGTHENED && jitterRealign.Type != SHORTENED)
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

        // add to jitter penalty as a function of the number of repeats found
        if(jitterRealign.Type == LENGTHENED || jitterRealign.Type == SHORTENED)
        {
            mJitterPenalty += jitterPenalty(qualityConfig, jitterRealign.RepeatCount);

            if(jitterRealign.Type == LENGTHENED)
                mLengthened++;
            else
                mShortened++;
        }

        return matchType;
    }

    private RawContext createRawContextFromCoreMatch(final SAMRecord record)
    {
        // check for an exact core match by which to centre the read
        if(mMaxCandidateDeleteLength < 5)
            return RawContext.INVALID_CONTEXT;

        int scLenLeft = record.getCigar().isLeftClipped() ? record.getCigar().getFirstCigarElement().getLength() : 0;
        int scLenRight = record.getCigar().isRightClipped() ? record.getCigar().getLastCigarElement().getLength() : 0;

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
                true, false, true, baseQuality, 0);
    }

    private boolean readHasIndelInCore(final SAMRecord record)
    {
        if(!record.getCigar().containsOperator(D) && !record.getCigar().containsOperator(I))
            return false;

        int variantLeftCorePos = mVariant.position() - (mReadContext.indexedBases().Index - mReadContext.indexedBases().LeftCoreIndex);
        int variantRightCorePos = mVariant.position() + (mReadContext.indexedBases().RightCoreIndex - mReadContext.indexedBases().Index);

        int currentPos = record.getAlignmentStart() - 1;
        int currentIndex = 0;

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
                currentIndex += element.getLength();
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
            IndexedBases readBases = new IndexedBases(position(), realignLeftReadIndex, record.getReadBases());

            ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                    readBases, record.getBaseQualities(), false, 0);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL)
                return new RealignedContext(EXACT, mReadContext.indexedBases().length(), realignLeftReadIndex);
        }

        int realignRightReadIndex = calcRightAlignmentIndex(record);

        if(realignRightReadIndex >= 0)
        {
            // still need to test even if this index matches the original readIndex since if the readIndex was in a delete
            // it will be have skipped above
            IndexedBases readBases = new IndexedBases(position(), realignRightReadIndex, record.getReadBases());

            ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                    readBases, record.getBaseQualities(), false, 0);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL)
                return new RealignedContext(RealignedType.EXACT, mReadContext.indexedBases().length(), realignRightReadIndex);
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
            final int[] existingCounts = mLpsCounts.get(index);
            if(readCount + allocCount > existingCounts[0] + existingCounts[0])
                break;

            ++index;
        }

        mLocalPhaseSets.add(index, lps);
        mLpsCounts.add(index, new int[] { readCount, (int)allocCount } );
    }

    private void countStrandedness(final SAMRecord record)
    {
        if(record.getFirstOfPairFlag())
            mForwardStrand++;
        else
            mReverseStrand++;
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

    private int qualityJitterPenalty() { return (int) mJitterPenalty; }

    private void checkImproperCount(final SAMRecord record)
    {
        if(!record.getReadPairedFlag() || !record.getProperPairFlag() || record.getSupplementaryAlignmentFlag())
        {
            mImproperPair++;
        }
    }
}
