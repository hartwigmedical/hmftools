package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.Collections;
import java.util.List;

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
    public final ChrBaseRegion RefLocation;
    public final int MapQual;
    public final int NMatches;
    public final int Score;
    public final int Flags;
    public final String Cigar;
    public final String XaTag;
    public final String MdTag;

    private final List<CigarElement> mCigarElements;
    private final int mSoftClipLeft;
    private final int mSoftClipRight;

    private final Orientation mOrientation;
    private final int mAlignedBases;

    private final int mRawSequenceStart;
    private final int mRawSequenceEnd;
    private int mSequenceStart;
    private int mSequenceEnd;
    private boolean mIsRequeried;

    private int mAdjustedAlignment;

    public AlignData(
            final ChrBaseRegion refLocation, final int sequenceStart, final int sequenceEnd, final int mapQual,
            final int score, final int flags, final String cigar, final int nMatches, final String xaTag, final String mdTag)
    {
        RefLocation = refLocation;
        mRawSequenceStart = sequenceStart;
        mRawSequenceEnd = sequenceEnd;
        MapQual = mapQual;
        NMatches = nMatches;
        Score = score;
        Flags = flags;

        Cigar = cigar;
        XaTag = xaTag;
        MdTag = mdTag;

        mCigarElements = CigarUtils.cigarElementsFromStr(cigar);

        mSoftClipLeft = mCigarElements.get(0).getOperator() == S ? mCigarElements.get(0).getLength() : 0;
        int lastIndex = mCigarElements.size() - 1;
        mSoftClipRight = mCigarElements.get(lastIndex).getOperator() == S ? mCigarElements.get(lastIndex).getLength() : 0;

        mOrientation = SamRecordUtils.isFlagSet(Flags, READ_REVERSE_STRAND) ? REVERSE : FORWARD;
        mAlignedBases = mCigarElements.stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        mSequenceStart = sequenceStart;
        mSequenceEnd = max(sequenceEnd - 1, 0);
        mIsRequeried = false;
        mAdjustedAlignment = mAlignedBases;
    }

    public Orientation orientation() { return mOrientation; }
    public boolean isForward() { return mOrientation.isForward(); }
    public boolean isReverse() { return mOrientation.isReverse(); }
    public boolean isSupplementary() { return SamRecordUtils.isFlagSet(Flags, SUPPLEMENTARY_ALIGNMENT); }

    public int maxSoftClipLength() { return max(mSoftClipLeft, mSoftClipRight); }
    public int leftSoftClipLength() { return mSoftClipLeft; }
    public int rightSoftClipLength() { return mSoftClipRight; }
    public int alignedBases() { return mAlignedBases; }
    public int segmentLength() { return mSequenceEnd - mSequenceStart + 1; }

    public void setFullSequenceData(final String fullSequence, final int fullSequenceLength)
    {
        if(mIsRequeried)
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

    public void setRequeriedSequenceCoords(int sequenceStart, int sequenceEnd)
    {
        mSequenceStart = sequenceStart;
        mSequenceEnd = sequenceEnd;
        mIsRequeried = true;
    }

    public boolean isRequeried() { return mIsRequeried; }

    public int sequenceStart() { return mSequenceStart; }
    public int sequenceEnd() { return mSequenceEnd; }

    public int rawSequenceStart() { return mRawSequenceStart; }
    public int rawSequenceEnd() { return mRawSequenceEnd; }

    public List<CigarElement> cigarElements() { return mCigarElements; }

    public int adjustedAlignment() { return mAdjustedAlignment; }

    public void setAdjustedAlignment(final String fullSequence, int inexactHomologyStart, int inexactHomologyEnd)
    {
        int sequenceLength = mSequenceEnd - mSequenceStart + 1;
        int scoreAdjustment = sequenceLength - Score;

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
    }

    public double calcModifiedMapQual()
    {
        double lengthFactor = pow(mAdjustedAlignment/ 100.0, 2);
        return (int)round(MapQual * min(1, lengthFactor));

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

    public List<AlternativeAlignment> altAlignments()
    {
        if(XaTag == null || XaTag.isEmpty())
            return Collections.emptyList();

        List<AlternativeAlignment> alternativeAlignments = AlternativeAlignment.fromLocationTag(XaTag);

        alternativeAlignments.add(
                0, new AlternativeAlignment(RefLocation.Chromosome, RefLocation.start(), orientation(), Cigar, MapQual));

        return alternativeAlignments;
    }

    public boolean isLowerAlignment(final AlignData other)
    {
        // returns true if this alignment is lower in genome position terms than the comparison alignment
        return RefLocation.compareTo(other.RefLocation) < 0;
    }

    public String toString()
    {
        return format("%s %s %s seq(%d-%d adj=%d-%d) score(%d) flags(%d) mapQual(%d align=%d adj=%d)",
                RefLocation, Cigar, mOrientation.isForward() ? "fwd" : "rev",  mRawSequenceStart, mRawSequenceEnd,
                mSequenceStart, mSequenceEnd, Score, Flags, MapQual, mAlignedBases, mAdjustedAlignment);
    }
}
