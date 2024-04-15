package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConfig.READ_ID_TRIMMER;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;

import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.READ_PAIRED;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
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

    private final int mJunctionReadIndex; // index within this read of the junction position
    private int mJunctionAssemblyIndex; // index within this read's junction assembly if the read's start position
    private int mLinkedAssemblyIndex; // index within this read's full linked assembly (if exists) if the read's start position

    // those past the junction
    private int mJunctionMatches;
    private int mJunctionMismatches;
    private int mReferenceMismatches;

    private Read mRead; // expect to be null unless required for BAM or read TSV writing

    public SupportRead(final Read read, final SupportType type, final int junctionReadIndex, final int matches, final int mismatches)
    {
        mType = type;

        mId = READ_ID_TRIMMER.trim(read.id());
        mChromosome = read.chromosome();
        mAlignmentStart = read.alignmentStart();
        mAlignmentEnd = read.alignmentEnd();
        mFlags = read.getFlags();
        mCigar = read.cigarString();

        mUnclippedStart = read.indelImpliedAlignmentStart() > 0 ?
                min(read.indelImpliedAlignmentStart(), read.unclippedStart()) : read.unclippedStart();

        mUnclippedEnd = read.indelImpliedAlignmentEnd() > 0 ?
                max(read.indelImpliedAlignmentEnd(), read.unclippedEnd()) : read.unclippedEnd();

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
        mNumOfEvents = read.numberOfEvents();

        mJunctionReadIndex = junctionReadIndex;
        mJunctionMatches = matches;
        mJunctionMismatches = mismatches;
        mReferenceMismatches = 0;
        mJunctionAssemblyIndex = -1;
        mLinkedAssemblyIndex = -1;

        mRead = read;
    }

    public SupportType type() { return mType; }

    public String id() { return mId; }
    public String chromosome() { return mChromosome; }
    public int alignmentStart() { return mAlignmentStart; }
    public int alignmentEnd() { return mAlignmentEnd; }
    public int unclippedStart() { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }
    public String mateChromosome() { return mMateChromosome; }
    public int mateAlignmentStart() { return mMateAlignmentStart; }
    public int mateAlignmentEnd() { return mMateAlignmentEnd; }
    public int baseLength() { return mBaseLength; }
    public int insertSize() { return mInsertSize; }
    public int trimCount() { return mTrimCount; }
    public String cigar() { return mCigar; }
    public byte orientation() { return isFlagSet(READ_REVERSE_STRAND) ? NEG_ORIENT : POS_ORIENT; }
    public SupplementaryReadData supplementaryData() { return mSupplementaryData; }
    public int mapQual() { return mMapQual; }
    public int numOfEvents() { return mNumOfEvents; }

    public boolean isReference() { return mIsReference; }
    public int sampleIndex() { return mSampleIndex; }

    public int flags() { return mFlags; }
    public boolean isSupplementary() { return isFlagSet(SUPPLEMENTARY_ALIGNMENT); }
    public boolean isPairedRead() { return isFlagSet(READ_PAIRED); }
    public boolean isMateUnmapped() { return isFlagSet(MATE_UNMAPPED); }
    public boolean isMateMapped() { return isFlagSet(READ_PAIRED) && !isFlagSet(MATE_UNMAPPED); }
    public boolean isDiscordant() { return mIsDiscordant; }

    public boolean isFlagSet(final SAMFlag flag) { return SamRecordUtils.isFlagSet(mFlags, flag); }

    public boolean matchesFragment(final SupportRead other) { return mId.equals(other.id()); }

    public void clearCachedRead() { mRead = null; }

    @Nullable
    public Read cachedRead() { return mRead; }

    public int junctionReadIndex() { return mJunctionReadIndex; }

    public void setJunctionAssemblyIndex(int index) { mJunctionAssemblyIndex = index; }
    public int junctionAssemblyIndex() { return mJunctionAssemblyIndex; }

    public void setLinkedAssemblyIndex(int index) { mLinkedAssemblyIndex = index; }
    public int linkedAssemblyIndex() { return mLinkedAssemblyIndex; }

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
        return format("type(%s) read(%s %d-%d %s) index(junc=%d asm=%d linked=%s) hqMatch(%d) mismatch(junc=%d ref=%d)",
                mType, mId, mAlignmentStart, mAlignmentEnd, mCigar, mJunctionReadIndex, mJunctionAssemblyIndex, mLinkedAssemblyIndex,
                mJunctionMatches, mJunctionMismatches, mReferenceMismatches);
    }
}
