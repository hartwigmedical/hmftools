package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvUtils.isDiscordant;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;

import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.READ_PAIRED;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.esvee.assembly.read.Read;

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

    // private final int mUnclippedStart;
    // private final int mUnclippedEnd;
    // private final int mNumberOfEvents;

    private final String mMateChromosome;
    private final int mMateAlignmentStart;
    private final int mMateAlignmentEnd;
    private final int mFlags;
    private final int mInsertSize;

    private final boolean mIsDiscordant;
    private final boolean mIsReference;
    private final int mSampleIndex;

    // fragment state
    private final SupplementaryReadData mSupplementaryData;
    private IndelCoords mIndelCoords;
    // private int mIndelImpliedAlignmentStart;
    // private int mIndelImpliedAlignmentEnd;
    private final int mMapQual;
    private final int mNumOfEvents;
    private final int mTrimCount;

    private final int mJunctionReadIndex; // index within this read of the junction position

    // those past the junction
    private int mJunctionMatches;
    private int mJunctionMismatches;
    private int mReferenceMismatches;

    private Read mRead; // expect to be null unless required for BAM or read TSV writing

    public SupportRead(final Read read, final SupportType type, final int junctionReadIndex, final int matches, final int mismatches)
    {
        mType = type;

        mId = read.id();
        mChromosome = read.chromosome();
        mAlignmentStart = read.alignmentStart();
        mAlignmentEnd = read.alignmentEnd();
        mFlags = read.getFlags();
        mCigar = read.cigarString();

        mMateChromosome = read.mateChromosome();
        mMateAlignmentStart = read.mateAlignmentStart();
        mMateAlignmentEnd = read.mateAlignmentEnd();
        mIsReference = read.isReference();
        mSampleIndex = read.sampleIndex();
        mIsDiscordant = isDiscordantFragment(read);
        mSupplementaryData = read.supplementaryData();
        mTrimCount = read.baseTrimCount();
        mInsertSize = read.insertSize();
        mMapQual = read.mappingQuality();
        mNumOfEvents = read.numberOfEvents();

        mJunctionReadIndex = junctionReadIndex;
        mJunctionMatches = matches;
        mJunctionMismatches = mismatches;
        mReferenceMismatches = 0;

        mRead = read;
    }

    public SupportType type() { return mType; }

    public String id() { return mId; }
    public String chromosome() { return mChromosome; }
    public int alignmentStart() { return mAlignmentStart; }
    public int alignmentEnd() { return mAlignmentEnd; }
    public String mateChromosome() { return mMateChromosome; }
    public int mateAlignmentStart() { return mMateAlignmentStart; }
    public int mateAlignmentEnd() { return mMateAlignmentEnd; }
    public String cigar() { return mCigar; }
    public boolean isReference() { return mIsReference; }
    public int sampleIndex() { return mSampleIndex; }
    public boolean isSupplementary() { return isFlagSet(SUPPLEMENTARY_ALIGNMENT); }
    public boolean isPairedRead() { return isFlagSet(READ_PAIRED); }
    public boolean isMateUnmapped() { return isFlagSet(MATE_UNMAPPED); }
    public boolean isDiscordant() { return mIsDiscordant; }

    public boolean isFlagSet(final SAMFlag flag) { return (mFlags & flag.intValue()) != 0; }
    public byte orientation() { return isFlagSet(READ_REVERSE_STRAND) ? NEG_ORIENT : POS_ORIENT; }
    public SupplementaryReadData supplementaryData() { return mSupplementaryData; }
    public int flags() { return mFlags; }
    public int insertSize() { return mTrimCount; }
    public int mapQual() { return mMapQual; }
    public int numOfEvents() { return mNumOfEvents; }
    public int trimCount() { return mTrimCount; }

    public boolean matchesFragment(final SupportRead other) { return mId.equals(other.id()); }

    public void clearCachedRead() { mRead = null; }

    @Nullable
    public Read cachedRead() { return mRead; }

    public int junctionReadIndex() { return mJunctionReadIndex; }

    public int mismatchCount() { return mJunctionMismatches + mReferenceMismatches; }
    public int junctionMismatches() { return mJunctionMismatches; }
    public int junctionMatches() { return mJunctionMatches; }
    public int referenceMismatches() { return mReferenceMismatches; }
    public void setReferenceMismatches(int mismatches) { mReferenceMismatches = mismatches; }

    public static boolean hasMatchingFragment(final List<SupportRead> support, final SupportRead read)
    {
        return support.stream().anyMatch(x -> x.matchesFragment(read));
    }

    public static List<SupportRead> findMatchingFragmentSupport(final List<SupportRead> support, final SupportRead read)
    {
        return support.stream().filter(x -> x.matchesFragment(read)).collect(Collectors.toList());
    }

    public String toString()
    {
        return format("type(%s) read(%s %d-%d %s) juncIndex(%d) hqMatch(%d) mismatch(junc=%d ref=%d)",
                mType, mId, mAlignmentStart, mAlignmentEnd, mCigar,
                mJunctionReadIndex, mJunctionMatches, mJunctionMismatches, mReferenceMismatches);
    }
}
