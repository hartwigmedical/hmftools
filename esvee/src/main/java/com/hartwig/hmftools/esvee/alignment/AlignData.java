package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

import htsjdk.samtools.CigarElement;

public class AlignData
{
    private final ChrBaseRegion mRefLocation;
    private final int mMapQual;
    private final int mNMatches;
    private final int mScore;
    private final int mFlags;
    private final String mCigar;
    private final String mXaTag;
    private final String mMdTag;

    private final List<CigarElement> mCigarElements;
    private final int mSoftClipLeft;
    private final int mSoftClipRight;

    private final Orientation mOrientation;
    private final int mAlignedBases;

    private final int mRawSequenceStart;
    private final int mRawSequenceEnd;
    private int mSequenceStart;
    private int mSequenceEnd;
    private boolean mDroppedOnRequery;
    private boolean mIsRequery;

    private int mAdjustedAlignment;
    private int mModifiedMapQual;

    private List<AlternativeAlignment> mAltAlignments; // lazy instantiation
    private boolean mUseLowMapQualAlignment;
    private AlternativeAlignment mLinkedAltAlignment;

    public AlignData(
            final ChrBaseRegion refLocation, final int sequenceStart, final int sequenceEnd, final int mapQual,
            final int score, final int flags, final String cigar, final int nMatches, final String xaTag, final String mdTag)
    {
        mRefLocation = refLocation;
        mRawSequenceStart = sequenceStart;
        mRawSequenceEnd = sequenceEnd;
        mMapQual = mapQual;
        mNMatches = nMatches;
        mScore = score;
        mFlags = flags;

        mCigar = cigar;
        mXaTag = xaTag == null ? "" : xaTag;
        mMdTag = mdTag;

        mCigarElements = CigarUtils.cigarElementsFromStr(cigar);

        mSoftClipLeft = mCigarElements.get(0).getOperator() == S ? mCigarElements.get(0).getLength() : 0;
        int lastIndex = mCigarElements.size() - 1;
        mSoftClipRight = mCigarElements.get(lastIndex).getOperator() == S ? mCigarElements.get(lastIndex).getLength() : 0;

        mOrientation = SamRecordUtils.isFlagSet(mFlags, READ_REVERSE_STRAND) ? REVERSE : FORWARD;
        mAlignedBases = mCigarElements.stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        mSequenceStart = sequenceStart;
        mSequenceEnd = max(sequenceEnd - 1, 0);
        mIsRequery = false;
        mDroppedOnRequery = false;
        mAdjustedAlignment = mAlignedBases;
        mModifiedMapQual = 0;

        mAltAlignments = null;
        mUseLowMapQualAlignment = false;
        mLinkedAltAlignment = null;
    }

    public ChrBaseRegion refLocation() { return mRefLocation; }
    public int mapQual() { return mMapQual; }
    public int nMatches() { return mNMatches; }
    public int score() { return mScore; }
    public int flags() { return mFlags; }
    public String cigar() { return mCigar; }
    public String xaTag() { return mXaTag; }
    public String mdTag() { return mMdTag; }

    public Orientation orientation() { return mOrientation; }
    public boolean isForward() { return mOrientation.isForward(); }
    public boolean isReverse() { return mOrientation.isReverse(); }
    public boolean isSupplementary() { return SamRecordUtils.isFlagSet(mFlags, SUPPLEMENTARY_ALIGNMENT); }

    public int leftSoftClipLength() { return mSoftClipLeft; }
    public int rightSoftClipLength() { return mSoftClipRight; }
    public int alignedBases() { return mAlignedBases; }
    public int segmentLength() { return mSequenceEnd - mSequenceStart + 1; }

    public void setFullSequenceData(final String fullSequence, final int fullSequenceLength)
    {
        if(mIsRequery)
            return;

        if(mOrientation.isReverse())
        {
            int newSequenceStart = (fullSequenceLength - 1) - (mRawSequenceEnd - 1);
            int newSequenceEnd = (fullSequenceLength - 1) - mSequenceStart;

            mSequenceStart = max(newSequenceStart, 0);
            mSequenceEnd = newSequenceEnd;
        }

        if(mSequenceStart < 0 || mSequenceStart > mSequenceEnd || mSequenceEnd > fullSequence.length())
        {
            SV_LOGGER.error("alignment({}) invalid subsequence request({}-{}) vs fullSequenceLength({})",
                    toString(), mSequenceStart, mSequenceEnd, fullSequenceLength);
        }
    }

    public void markDroppedOnRequery() { mDroppedOnRequery = true; }
    public boolean droppedOnRequery() { return mDroppedOnRequery; }

    public void setRequeriedSequenceCoords(int sequenceStart, int sequenceEnd)
    {
        mSequenceStart = sequenceStart;
        mSequenceEnd = sequenceEnd;
        mIsRequery = true;
    }

    public boolean isRequeried() { return mIsRequery; }

    public int sequenceStart() { return mSequenceStart; }
    public int sequenceEnd() { return mSequenceEnd; }

    public int rawSequenceStart() { return mRawSequenceStart; }
    public int rawSequenceEnd() { return mRawSequenceEnd; }

    public List<CigarElement> cigarElements() { return mCigarElements; }

    public int adjustedAlignment() { return mAdjustedAlignment; }
    public double modifiedMapQual() { return mModifiedMapQual; }
    public boolean exceedsMapQualThreshold() { return mModifiedMapQual >= ALIGNMENT_MIN_MOD_MAP_QUAL; }

    public void setAdjustedAlignment(final String fullSequence, int inexactHomologyStart, int inexactHomologyEnd)
    {
        int sequenceLength = mSequenceEnd - mSequenceStart + 1;
        int scoreAdjustment = sequenceLength - mScore;

        int seqStart = inexactHomologyStart > 0 ? mSequenceStart + inexactHomologyStart : mSequenceStart;
        int seqEnd = inexactHomologyEnd > 0 ? mSequenceEnd - inexactHomologyEnd : mSequenceEnd;

        if(seqStart < 0 || seqStart > seqEnd || seqEnd >= fullSequence.length())
        {
            mAdjustedAlignment = 0;
            return;
        }

        String alignedBases = fullSequence.substring(seqStart, seqEnd + 1);
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(alignedBases.getBytes());

        int trimmedAlignmentLength = calcTrimmedBaseLength(0, alignedBases.length() - 1, repeats);
        mAdjustedAlignment = max(trimmedAlignmentLength - scoreAdjustment, 0);

        double lengthFactor = pow(mAdjustedAlignment/ (double)max(100, mAlignedBases), 2);
        mModifiedMapQual = (int)round(mMapQual * min(1, lengthFactor));
    }

    public void markUseLowMapQualAlignment() { mUseLowMapQualAlignment = true; }
    public boolean useLowMapQualAlignment() { return mUseLowMapQualAlignment; }

    public boolean hasLinkedAltAlignment() { return mLinkedAltAlignment != null; }
    public AlternativeAlignment linkedAltAlignment() { return mLinkedAltAlignment; }
    public void setLinkedAltAlignment(final AlternativeAlignment alignment) { mLinkedAltAlignment = alignment; }

    public static AlignData from(final BwaMemAlignment alignment, final RefGenomeVersion refGenomeVersion)
    {
        int chrIndex = alignment.getRefId();

        if(chrIndex < 0 || chrIndex >= HumanChromosome.values().length)
            return null;

        String chromosome = refGenomeVersion.versionedChromosome(HumanChromosome.values()[chrIndex].toString());

        // note the +1 to the ref start position
        return new AlignData(
                new ChrBaseRegion(chromosome, alignment.getRefStart() + 1, alignment.getRefEnd()),
                alignment.getSeqStart(), alignment.getSeqEnd(), alignment.getMapQual(), alignment.getAlignerScore(),
                alignment.getSamFlag(), alignment.getCigar(), alignment.getNMismatches(), alignment.getXATag(), alignment.getMDTag());
    }

    public boolean hasAltAlignments() { return !mXaTag.isEmpty(); }

    public List<AlternativeAlignment> altAlignments()
    {
        if(mXaTag.isEmpty())
            return Collections.emptyList();

        if(mAltAlignments != null)
            return mAltAlignments;

        mAltAlignments = AlternativeAlignment.fromLocationTag(mXaTag);

        return mAltAlignments;
    }

    public List<AlternativeAlignment> allAlignments()
    {
        List<AlternativeAlignment> altAlignments = altAlignments();

        if(altAlignments.isEmpty())
            return altAlignments;

        // otherwise include the default alignment as well
        List<AlternativeAlignment> allAlignments = Lists.newArrayListWithCapacity(altAlignments.size() + 1);

        allAlignments.add(new AlternativeAlignment(mRefLocation.Chromosome, mRefLocation.start(), orientation(), mCigar, mMapQual));
        allAlignments.addAll(altAlignments);
        return allAlignments;
    }

    public boolean isLowerAlignment(final AlignData other)
    {
        // returns true if this alignment is lower in genome position terms than the comparison alignment
        return mRefLocation.compareTo(other.mRefLocation) < 0;
    }

    public String toString()
    {
        return format("%s %s %s seq(%d-%d adj=%d-%d) score(%d) flags(%d) mapQual(%d adj=%d) aligned(%d adj=%d)",
                mRefLocation, mCigar, mOrientation.isForward() ? "fwd" : "rev",  mRawSequenceStart, mRawSequenceEnd,
                mSequenceStart, mSequenceEnd, mScore, mFlags, mMapQual, mModifiedMapQual, mAlignedBases, mAdjustedAlignment);
    }
}
