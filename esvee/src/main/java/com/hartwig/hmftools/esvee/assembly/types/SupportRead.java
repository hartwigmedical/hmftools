package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static htsjdk.samtools.SAMFlag.MATE_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.READ_PAIRED;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.READ_UNMAPPED;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.List;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFlag;

public class SupportRead
{
    private final SupportType mType;

    // inherited read properties & state
    private final String mId;

    private final String mChromosome;
    private final int mAlignmentStart;
    private final int mAlignmentEnd;
    private final String mCigar;

    private final int mUnclippedStart;
    private final int mUnclippedEnd;

    private final String mMateChromosome;
    private final int mMateAlignmentStart;
    private final int mMateAlignmentEnd;
    private final int mFlags;
    private final int mBaseLength;

    private final boolean mIsDiscordant;
    private final boolean mIsReference;
    private final int mSampleIndex;

    // fragment state
    private final SupplementaryReadData mSupplementaryData;
    private final int mMapQual;
    private final int mNumOfEvents;
    private final int mInsertSize;
    private final int mTrimCount;
    private final boolean mHasIndel;
    private IndelCoords mIndelCoords;
    private final boolean mHasLineTail;

    // the distance from the read's start (ie index not position) to the assembly junction index
    // if the read start is before the junction index then the value is negative
    private final int mJunctionReadStartDistance;

    private int mFullAssemblyIndex; // index within this read's full linked assembly sequence (if exists) if the read's start position
    private Orientation mFullAssemblyOrientation;
    private int mInferredFragmentLength;
    private SupportType mBreakendType;

    // those past the junction
    private int mJunctionMatches;
    private int mJunctionMismatches;
    private Integer mReferenceMismatches;

    private Read mRead; // expect to be null unless required for BAM or read TSV writing

    public SupportRead(final Read read, final SupportType type, final int junctReadStartDistance, final int matches, final int mismatches)
    {
        mType = type;

        mId = read.id();

        mChromosome = read.chromosome();
        mAlignmentStart = read.alignmentStart();
        mAlignmentEnd = read.alignmentEnd();
        mFlags = read.getFlags();
        mCigar = read.cigarString();

        mUnclippedStart = read.unclippedStart();
        mUnclippedEnd = read.unclippedEnd();

        mMateChromosome = read.mateChromosome();
        mMateAlignmentStart = read.mateAlignmentStart();
        mMateAlignmentEnd = read.mateAlignmentEnd();
        mIsReference = read.isReference();
        mSampleIndex = read.sampleIndex();
        mIsDiscordant = isDiscordantFragment(read);
        mSupplementaryData = read.supplementaryData();
        mBaseLength = read.basesLength();
        mInsertSize = abs(read.bamRecord().getInferredInsertSize());
        mTrimCount = read.baseTrimCount();
        mMapQual = read.mappingQuality();
        mNumOfEvents = read.numOfEvents();
        mHasIndel = read.indelCoords() != null;
        mIndelCoords = read.indelCoords() != null && read.indelCoords().Length >= MIN_INDEL_LENGTH ? read.indelCoords() : null;
        mHasLineTail = read.hasLineTail();

        mJunctionMatches = matches;
        mJunctionMismatches = mismatches;
        mReferenceMismatches =  null;

        mJunctionReadStartDistance = junctReadStartDistance;
        mFullAssemblyIndex = -1;
        mFullAssemblyOrientation = null;
        mInferredFragmentLength = -1;
        mBreakendType = null;

        mRead = read;
    }

    public SupportType type() { return mType; }

    public String id() { return mId; }
    public String chromosome() { return mChromosome; }
    public int alignmentStart() { return mAlignmentStart; }
    public int alignmentEnd() { return mAlignmentEnd; }
    public int unclippedStart() { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }
    public boolean isLeftClipped() { return mUnclippedStart < mAlignmentStart; }
    public boolean isRightClipped() { return mUnclippedEnd > mAlignmentEnd; }
    public int leftClipLength() { return max(mAlignmentStart - mUnclippedStart, 0); }
    public int rightClipLength() { return max(mUnclippedEnd - mAlignmentEnd, 0); }
    public String mateChromosome() { return mMateChromosome; }
    public int mateAlignmentStart() { return mMateAlignmentStart; }
    public int mateAlignmentEnd() { return mMateAlignmentEnd; }
    public int baseLength() { return mBaseLength; }
    public int insertSize() { return mInsertSize; }
    public int trimCount() { return mTrimCount; }
    public boolean hasIndel() { return mHasIndel; }
    public IndelCoords indelCoords() { return mIndelCoords; }
    public String cigar() { return mCigar; }
    public Orientation orientation() { return isFlagSet(READ_REVERSE_STRAND) ? REVERSE : FORWARD; }
    public Orientation mateOrientation() { return isFlagSet(MATE_REVERSE_STRAND) ? REVERSE : FORWARD; }
    public SupplementaryReadData supplementaryData() { return mSupplementaryData; }
    public int mapQual() { return mMapQual; }
    public int numOfEvents() { return mNumOfEvents; }
    public boolean hasLineTail() { return mHasLineTail; }

    public boolean isReference() { return mIsReference; }
    public int sampleIndex() { return mSampleIndex; }

    public int flags() { return mFlags; }
    public boolean isSupplementary() { return isFlagSet(SUPPLEMENTARY_ALIGNMENT); }
    public boolean isPairedRead() { return isFlagSet(READ_PAIRED); }
    public boolean firstInPair() { return isFlagSet(FIRST_OF_PAIR); }
    public boolean isUnmapped() { return isFlagSet(READ_UNMAPPED); }
    public boolean isMateUnmapped() { return isFlagSet(MATE_UNMAPPED); }
    public boolean isMateMapped() { return isFlagSet(READ_PAIRED) && !isFlagSet(MATE_UNMAPPED); }

    public boolean isDiscordant() { return mIsDiscordant; }

    public boolean isFlagSet(final SAMFlag flag) { return SamRecordUtils.isFlagSet(mFlags, flag); }

    public int junctionMismatches() { return mJunctionMismatches; }
    public int junctionMatches() { return mJunctionMatches; }

    public int referenceMismatches() { return mReferenceMismatches != null ? mReferenceMismatches : -1; }
    public boolean hasReferenceMismatches() { return mReferenceMismatches != null; }

    public void setReferenceMismatches(int mismatches) { mReferenceMismatches = mismatches; }

    @Nullable
    public Read cachedRead() { return mRead; }

    public void clearCachedRead() { mRead = null; }

    public int junctionReadStartDistance() { return mJunctionReadStartDistance; }

    public void setFullAssemblyInfo(int assemblyIndex, final Orientation orientation)
    {
        mFullAssemblyIndex = assemblyIndex;
        mFullAssemblyOrientation = orientation;
    }

    public int fullAssemblyIndexStart() { return mFullAssemblyIndex; }
    public int fullAssemblyIndexEnd() { return mFullAssemblyIndex + mBaseLength - 1; }
    public Orientation fullAssemblyOrientation() { return mFullAssemblyOrientation; }

    public int inferredFragmentLength() { return mInferredFragmentLength; }
    public void setInferredFragmentLength(int length) { mInferredFragmentLength = length; }

    public SupportType breakendSupportType() { return mBreakendType; }
    public void setBreakendSupportType(final SupportType breakendType) { mBreakendType = breakendType; }

    public boolean matchesFragment(final SupportRead other, boolean allowReadMatch)
    {
        if(!mId.equals(other.id()))
            return false;

        return allowReadMatch || mFlags != other.flags();
    }

    public boolean matchesFragment(final Read other, boolean allowReadMatch)
    {
        if(!mId.equals(other.id()))
            return false;

        return allowReadMatch || mFlags != other.getFlags();
    }

    public static boolean hasFragmentOtherRead(final List<SupportRead> support, final SupportRead read)
    {
        return hasFragmentOtherRead(support, read.cachedRead());
    }

    public static boolean hasFragmentOtherRead(final List<SupportRead> support, final Read read)
    {
        return support.stream().anyMatch(x -> x.matchesFragment(read, false));
    }

    public static boolean hasMatchingFragmentRead(final List<SupportRead> support, final SupportRead read)
    {
        return support.stream().anyMatch(x -> x.matchesFragment(read, true));
    }

    public static boolean hasMatchingFragmentRead(final List<SupportRead> support, final Read read)
    {
        return support.stream().anyMatch(x -> x.matchesFragment(read, true));
    }

    public String toString()
    {
        return format("type(%s) read(%s %s:%d-%d %s %d) index(juncDist=%d asm=%d:%d) hqMatch(%d) mismatch(junc=%d ref=%d)",
                mType, mId, mChromosome, mAlignmentStart, mAlignmentEnd, mCigar, orientation().asByte(), mJunctionReadStartDistance,
                mFullAssemblyIndex, mFullAssemblyOrientation != null ? mFullAssemblyOrientation.asByte() : 0,
                mJunctionMatches, mJunctionMismatches, mReferenceMismatches != null ? mReferenceMismatches : -1);
    }
}
