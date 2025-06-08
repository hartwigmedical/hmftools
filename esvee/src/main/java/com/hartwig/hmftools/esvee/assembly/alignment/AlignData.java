package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_CALC_SCORE_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_CALC_SCORE_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_INDEL_RECOVERY_MIN_MAP_QUAL;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.common.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
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
    private ChrBaseRegion mRefLocation;
    private final int mMapQual;
    private final int mNMatches;
    private final int mScore;
    private final int mFlags;
    private final String mCigar;
    private final String mXaTag;
    private final String mMdTag;

    private final List<CigarElement> mCigarElements;
    private int mSoftClipLeft;
    private int mSoftClipRight;

    private Orientation mOrientation;
    private final int mAlignedBases;

    private int mRawSequenceStart;
    private int mRawSequenceEnd;
    private int mSequenceStart;
    private int mSequenceEnd;
    private boolean mDroppedOnRequery;
    private boolean mHasSequenceCoordsSet;

    private int mAdjustedAlignment;
    private int mModifiedMapQual;

    private List<AlternativeAlignment> mRawAltAlignments; // lazy instantiation of BWA-returned alt alignments
    private AlternativeAlignment mSelectedAltAlignment; // set if the default is not selected
    private List<AlternativeAlignment> mUnselectedAltAlignments; // those other than the one selected
    private boolean mHasLowMapQualShortSvLink;

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

        if(!mCigarElements.isEmpty())
        {
            mSoftClipLeft = mCigarElements.get(0).getOperator() == S ? mCigarElements.get(0).getLength() : 0;
            int lastIndex = mCigarElements.size() - 1;
            mSoftClipRight = mCigarElements.get(lastIndex).getOperator() == S ? mCigarElements.get(lastIndex).getLength() : 0;
        }
        else
        {
            // TEMP: should not occur
            mSoftClipLeft = mSoftClipRight = 0;
        }

        mOrientation = SamRecordUtils.isFlagSet(mFlags, READ_REVERSE_STRAND) ? REVERSE : FORWARD;
        mAlignedBases = mCigarElements.stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        mSequenceStart = sequenceStart;
        mSequenceEnd = max(sequenceEnd - 1, 0);
        mHasSequenceCoordsSet = false;
        mDroppedOnRequery = false;
        mAdjustedAlignment = mAlignedBases;
        mModifiedMapQual = 0;

        mRawAltAlignments = null;
        mSelectedAltAlignment = null;
        mUnselectedAltAlignments = null;
        mHasLowMapQualShortSvLink = false;
    }

    public ChrBaseRegion refLocation() { return mRefLocation; }
    public String chromosome() { return mRefLocation.Chromosome; }
    public int positionStart() { return mRefLocation.start(); }
    public int positionEnd() { return mRefLocation.end(); }
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

    public void setSoftClipLengths(int left, int right)
    {
        mSoftClipLeft = left;
        mSoftClipRight = right;
    }

    public void setFullSequenceData(final String fullSequence, final int fullSequenceLength)
    {
        if(mHasSequenceCoordsSet)
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

        mHasSequenceCoordsSet = true;
    }

    public void markDroppedOnRequery() { mDroppedOnRequery = true; }
    public boolean droppedOnRequery() { return mDroppedOnRequery; }

    public void setRequeriedSequenceCoords(int sequenceStart, int sequenceEnd)
    {
        mSequenceStart = sequenceStart;
        mSequenceEnd = sequenceEnd;
        mHasSequenceCoordsSet = true;
    }

    public int sequenceStart() { return mSequenceStart; }
    public int sequenceEnd() { return mSequenceEnd; }

    public int rawSequenceStart() { return mRawSequenceStart; }
    public int rawSequenceEnd() { return mRawSequenceEnd; }

    public List<CigarElement> cigarElements() { return mCigarElements; }

    public int adjustedAlignment() { return mAdjustedAlignment; }
    public double modifiedMapQual() { return mModifiedMapQual; }

    private int calcAdjustedIndelScore()
    {
        if(mMapQual < ALIGNMENT_INDEL_RECOVERY_MIN_MAP_QUAL)
            return 0;

        int indelLength = 0;

        for(CigarElement element : mCigarElements)
        {
            if(element.getOperator() == S) // must be fully aligned reads, ie no soft-clips
                return 0;

            if(element.getOperator().isIndel() && element.getLength() >= MIN_VARIANT_LENGTH)
            {
                if(indelLength != 0) // only allow 1 long INDEL, so exit if a second is found
                    return 0;

                indelLength = element.getLength();

                if(element.getOperator() == D)
                    indelLength = -indelLength;
            }
        }

        return indelLength;
    }

    public void setAdjustedAlignment(final String fullSequence, int inexactHomologyStart, int inexactHomologyEnd)
    {
        int sequenceLength = mSequenceEnd - mSequenceStart + 1;

        // adjust the score for CIGAR indels >= min length from single alignments
        int adjustedScore = mScore;

        int adjustedIndelScore = calcAdjustedIndelScore();

        if(adjustedIndelScore != 0)
        {
            if(adjustedIndelScore < 0)
                sequenceLength += adjustedIndelScore; // remove deleted segment
            else
                sequenceLength -= adjustedIndelScore;

            adjustedScore += abs(adjustedIndelScore) + BWA_GAP_OPEN_PENALTY;
        }

        int scoreAdjustment = sequenceLength - adjustedScore;

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

        int inexactHomologyLength = inexactHomologyStart + inexactHomologyEnd;
        double adjustedSequenceLength = max(sequenceLength - inexactHomologyLength, 0);

        if((adjustedScore + ALIGNMENT_CALC_SCORE_FACTOR) / adjustedSequenceLength < ALIGNMENT_CALC_SCORE_THRESHOLD)
        {
            mModifiedMapQual = 0;
        }
        else
        {
            double lengthFactor = pow(mAdjustedAlignment / (double) max(100, mAlignedBases - inexactHomologyLength), 2);
            mModifiedMapQual = (int) round(mMapQual * min(1, lengthFactor));
        }
    }

    public boolean hasAltAlignments() { return !mXaTag.isEmpty(); }

    public List<AlternativeAlignment> rawAltAlignments()
    {
        if(mXaTag.isEmpty())
            return Collections.emptyList();

        if(mRawAltAlignments != null)
            return mRawAltAlignments;

        mRawAltAlignments = AlternativeAlignment.fromLocationTag(mXaTag);

        return mRawAltAlignments;
    }

    public List<AlternativeAlignment> allAlignments()
    {
        List<AlternativeAlignment> altAlignments = rawAltAlignments();

        AlternativeAlignment defaultAlignment = new AlternativeAlignment(
                mRefLocation.Chromosome, mRefLocation.start(), orientation(), mCigar, mMapQual);

        if(altAlignments.isEmpty())
        {
            if(mRefLocation.Chromosome.startsWith(CHR_PREFIX))
            {
                if(MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38.stream().noneMatch(x -> x.overlaps(mRefLocation)))
                    return altAlignments;
            }
            else
            {
                if(MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37.stream().noneMatch(x -> x.overlaps(mRefLocation)))
                    return altAlignments;
            }
        }

        // otherwise include the default alignment as well
        List<AlternativeAlignment> allAlignments = Lists.newArrayListWithCapacity(altAlignments.size() + 1);

        allAlignments.add(defaultAlignment);
        allAlignments.addAll(altAlignments);
        return allAlignments;
    }

    public boolean hasLowMapQualAlignment() { return mSelectedAltAlignment != null || mUnselectedAltAlignments != null; }
    public boolean hasSelectedAltAlignment() { return mSelectedAltAlignment != null; }
    public boolean hasLowMapQualShortSvLink() { return mHasLowMapQualShortSvLink; }
    public AlternativeAlignment selectedAltAlignment() { return mSelectedAltAlignment; }

    public void setSelectedAltAlignments(
            final AlternativeAlignment selectedAlignment, final List<AlternativeAlignment> otherAlignments, boolean hasShortSvLink)
    {
        mSelectedAltAlignment = selectedAlignment;
        mUnselectedAltAlignments = otherAlignments;
        mHasLowMapQualShortSvLink = hasShortSvLink;
    }

    public void setSelectedLowMapQualAltAlignment(final AlternativeAlignment selectedAlignment)
    {
        mSelectedAltAlignment = selectedAlignment;
        mModifiedMapQual = ALIGNMENT_MIN_MOD_MAP_QUAL;
        mRawAltAlignments = null;
        mSoftClipRight = mSoftClipLeft = 0;
    }

    public List<AlternativeAlignment> unselectedAltAlignments()
    {
        return mUnselectedAltAlignments != null ? mUnselectedAltAlignments : Collections.emptyList();
    }

    public boolean isLowerAlignment(final AlignData other)
    {
        // returns true if this alignment is lower in genome position terms than the comparison alignment
        return mRefLocation.compareTo(other.mRefLocation) < 0;
    }

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

    public String toString()
    {
        return format("%s:%d %s seq(%d-%d adj=%d-%d) score(%d) flags(%d) mapQual(%d adj=%d) aligned(%d adj=%d) md(%s)",
                mRefLocation, mOrientation.asByte(), mCigar, mRawSequenceStart, mRawSequenceEnd,
                mSequenceStart, mSequenceEnd, mScore, mFlags, mMapQual, mModifiedMapQual, mAlignedBases, mAdjustedAlignment, mMdTag);
    }
}
