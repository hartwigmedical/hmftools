package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.esvee.alignment.AlignmentOutcome.NO_SET;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.calcTrimmedRefBaseLength;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.readQualFromJunction;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.convertedIndelCrossesJunction;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findInsertedBases;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.REMOTE_REGION;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;

import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.alignment.AlignmentOutcome;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public class JunctionAssembly
{
    private int mAssemblyId;
    private final Junction mJunction;

    private int mJunctionIndex; // position of the junction in the read bases

    // aligned position on the ref base side
    private int mRefBasePosition;
    private final List<RefBaseIndel> mRefBaseIndels;

    private final IndelCoords mIndelCoords;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<SupportRead> mSupport;
    private final List<Read> mCandidateSupport;
    private final List<Read> mUnmappedCandidates;
    private final List<RefSideSoftClip> mRefSideSoftClips;

    private final List<RepeatInfo> mRepeatInfo;
    private String mRefBasesRepeatedTrimmed;
    private int mRefBaseTrimLength;

    private final List<RemoteRegion> mRemoteRegions;

    private PhaseGroup mPhaseGroup;
    private AssemblyOutcome mOutcome;
    private AlignmentOutcome mAlignmentOutcome;
    private String mAssemblyAlignmentInfo;

    // info only
    private final String mInitialReadId;
    private int mMergedAssemblies;
    private int mMismatchReadCount;

    private final AssemblyStats mStats;

    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] baseQualities, final List<SupportRead> assemblySupport,
            final List<RepeatInfo> repeatInfo)
    {
        mJunction = junction;

        // set initial bounds from the assembly support
        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        IndelCoords indelCoords = null;
        mRefBasePosition = junction.Position; // initialised to the same prior to extending ref bases

        for(SupportRead support : assemblySupport)
        {
            Read read = support.cachedRead();

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }

            if(read.indelCoords() != null && indelCoords == null)
            {
                indelCoords = read.indelCoords();

                if(indelCoords.isInsert())
                    indelCoords.setInsertedBases(findInsertedBases(read));
            }
        }

        mIndelCoords = indelCoords;

        mInitialReadId = maxJunctionBaseQualRead != null ? maxJunctionBaseQualRead.id() :
                (!assemblySupport.isEmpty() ? assemblySupport.get(0).id() : "null");

        mBases = bases;
        mBaseQuals = baseQualities;

        mJunctionIndex = junction.isForward() ? 0 : mBases.length - 1;

        mRefBaseIndels = Lists.newArrayList();
        mSupport = Lists.newArrayList(assemblySupport);
        mCandidateSupport = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();
        mRepeatInfo = repeatInfo;
        mRefSideSoftClips = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mMergedAssemblies = 0;
        mPhaseGroup = null;
        mOutcome = UNSET;
        mAlignmentOutcome = NO_SET;
        mAssemblyAlignmentInfo = null;
        mMismatchReadCount = 0;

        mStats = new AssemblyStats();
        assemblySupport.forEach(x -> mStats.addRead(x, mJunction, x.cachedRead()));
    }

    public void setId(int id) { mAssemblyId = id; }
    public int id() { return mAssemblyId; }

    public Junction junction() { return mJunction; }
    public boolean isForwardJunction() { return mJunction.isForward(); }

    public boolean indel() { return mJunction.IndelBased; }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionIndex; };
    public void setJunctionIndex(int index) { mJunctionIndex = index; };

    // eg 21 bases, junction index at 10 (so 0-9 = 10 before, 11-20 = 10 after), note: doesn't count the junction base
    public int lowerDistanceFromJunction() { return mJunctionIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionIndex - 1; };

    public int refBaseLength() { return (mJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction()) + 1; }

    // base count beyond the junction
    public int extensionLength() { return mJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    public int refBasePosition() { return mRefBasePosition; }
    public int baseLength() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public String initialReadId() { return mInitialReadId; }
    public IndelCoords indelCoords() { return mIndelCoords; }

    public List<SupportRead> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public AssemblyStats stats() { return mStats; }

    public int mismatchReadCount() { return mMismatchReadCount; }
    public void addMismatchReadCount(int count) { mMismatchReadCount += count; }

    public void addRead(
            final Read read, final ReadAssemblyIndices readAssemblyIndices, final SupportType type, @Nullable final SupportRead existingSupport)
    {
        if(readAssemblyIndices == ReadAssemblyIndices.INVALID_INDICES)
            return;

        int mismatchCount = 0;
        int highQualMatchCount = 0;
        int assemblyIndex = readAssemblyIndices.AssemblyIndexStart;

        for(int i = readAssemblyIndices.ReadIndexStart; i <= readAssemblyIndices.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length || i >= read.getBases().length)
            {
                // can possibly happen with INDELs in the CIGAR or not now that aligned positions are used to establish read coordinates?
                // SV_LOGGER.debug("i({}) readCoords({}) assembly({}) read({})", i, readCoords, toString(), read.toString());
                break;
            }

            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;

                if(qual >= LOW_BASE_QUAL_THRESHOLD)
                    ++highQualMatchCount;
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    if(qual >= LOW_BASE_QUAL_THRESHOLD)
                        ++highQualMatchCount;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                }
                else
                {
                    ++mismatchCount;
                }
            }
        }

        if(existingSupport == null)
        {
            addSupport(read, type, readAssemblyIndices.JunctionIndex, highQualMatchCount, mismatchCount, 0);
        }
        else
        {
            existingSupport.setReferenceMismatches(mismatchCount);
        }
    }

    public void addSupport(
            final Read read, final SupportType type, int readJunctionIndex, int matches, int mismatches, int refMismatches)
    {
        boolean isIndelCrossingJunction = convertedIndelCrossesJunction(mJunction, read);
        SupportType adjustedType = type == JUNCTION && isIndelCrossingJunction ? INDEL : type;
        SupportRead support = new SupportRead(read, adjustedType, readJunctionIndex, matches, mismatches);
        support.setReferenceMismatches(refMismatches);

        mSupport.add(support);
        mStats.addRead(support, mJunction, read);
    }

    public void extendRefBases(int newRefBasePosition, final List<RefBaseIndel> refBaseIndels, final RefGenomeInterface refGenome)
    {
        // extend the number of ref bases to accommodate new ref base information from existing or new reads
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        // TODO: factor in indels to base length and any ref genome bases
        // int totalIndelLength = refBaseIndels.stream().mapToInt(x -> x.Length).sum();
        int totalIndelLength = 0;

        int refBaseExtension = abs(newRefBasePosition - mRefBasePosition) + totalIndelLength;

        int newBaseLength = mBases.length + refBaseExtension;

        int existingRefBasePosition = mRefBasePosition;

        mRefBasePosition = newRefBasePosition;

        boolean isForwardJunction = mJunction.isForward();

        if(isForwardJunction)
            mJunctionIndex += refBaseExtension;

        int baseOffset = isForwardJunction ? refBaseExtension : 0;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        byte[] refGenomeBases = null;

        if(refGenome != null)
        {
            int refBaseStart, refBaseEnd;

            if(mJunction.isForward())
            {
                refBaseStart = newRefBasePosition;
                refBaseEnd = existingRefBasePosition - 1;
            }
            else
            {
                refBaseStart = existingRefBasePosition + 1;
                refBaseEnd = newRefBasePosition;
            }

            refGenomeBases = refGenome.getBases(mJunction.Chromosome, refBaseStart, refBaseEnd);
        }

        // forward junction: 0 up to base offset - 1 -> set to zero or ref genome bases, base offset to new length -> copy bases
        // reverse junction: 0 up to base offset - 1 -> copy bases, base offset to new length -> set to zero or ref genome bases

        int refBaseIndex = 0;

        for(int i = 0; i < newBaseLength; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    if(refGenomeBases != null && refBaseIndex < refGenomeBases.length)
                        mBases[i] = refGenomeBases[refBaseIndex];
                    else
                        mBases[i] = 0;

                    ++refBaseIndex;
                    mBaseQuals[i] = 0;
                }
                else
                {
                    mBases[i] = existingBases[i - baseOffset];
                    mBaseQuals[i] = existingQuals[i - baseOffset];
                }
            }
            else
            {
                if(i < existingBases.length)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                }
                else
                {
                    if(refGenomeBases != null && refBaseIndex < refGenomeBases.length)
                        mBases[i] = refGenomeBases[refBaseIndex];
                    else
                        mBases[i] = 0;

                    ++refBaseIndex;
                    mBaseQuals[i] = 0;
                }
            }
        }
    }

    public void expandExtensionBases(final byte[] extensionBases, final byte[] extensionBaseQuals, final List<SupportRead> supportReads)
    {
        // extend the same of ref bases to accommodate new ref base information from existing or new reads
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int newExtensionLength = extensionBases.length - 1; // since the extension bases include the junction index
        int baseLengthChange = newExtensionLength - extensionLength();
        int existingBaseLength = mBases.length;
        int newBaseLength = existingBaseLength + baseLengthChange;

        if(mJunction.isReverse())
            mJunctionIndex += baseLengthChange;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int newExtBaseIndex = 0;
        for(int i = 0; i < newBaseLength; ++i)
        {
            if(mJunction.isForward())
            {
                if(i < existingBaseLength)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                }
                else
                {
                    mBases[i] = extensionBases[newExtBaseIndex];
                    mBaseQuals[i] = extensionBaseQuals[newExtBaseIndex];
                    ++newExtBaseIndex;
                }
            }
            else
            {
                if(i < baseLengthChange)
                {
                    mBases[i] = extensionBases[newExtBaseIndex];
                    mBaseQuals[i] = extensionBaseQuals[newExtBaseIndex];
                    ++newExtBaseIndex;
                }
                else
                {
                    mBases[i] = existingBases[i - baseLengthChange];
                    mBaseQuals[i] = existingQuals[i - baseLengthChange];
                }
            }
        }

        for(SupportRead support : supportReads)
        {
            mSupport.add(support);
            mStats.addRead(support, mJunction, support.cachedRead());
        }
    }

    public void removeSupportReads(final Set<String> readIds)
    {
        int index = 0;

        while(index < mSupport.size())
        {
            if(readIds.contains(mSupport.get(index).id()))
                mSupport.remove(index);
            else
                ++index;
        }
    }

    public void trimRefBases(int newRefPosition)
    {
        // trims to the specified ref position, and assumes no indels are included
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int existingRefPosition = refBasePosition();

        int newBaseLength = mBases.length - abs(existingRefPosition - newRefPosition);
        int baseLengthChange = mBases.length - newBaseLength;

        if(mJunction.isForward())
            mJunctionIndex -= baseLengthChange;

        mRefBasePosition = newRefPosition;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int baseOffset = mJunction.isForward() ? baseLengthChange : 0;

        for(int i = 0; i < newBaseLength; ++i)
        {
            mBases[i] = existingBases[i + baseOffset];
            mBaseQuals[i] = existingQuals[i + baseOffset];
        }

        buildRepeatInfo();
    }

    public boolean hasReadSupport(final Read read)
    {
        return read != null && mSupport.stream().anyMatch(x -> x.cachedRead() == read);
    }

    public void clearSupportCachedReads() { mSupport.forEach(x -> x.clearCachedRead()); }

    // caching repeat info needs careful consideration since any extension of ref bases invalidates the values,
    // at least for +ve orientation assemblies
    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
        {
            mRepeatInfo.addAll(repeats);
            mRefBasesRepeatedTrimmed = RepeatInfo.buildTrimmedRefBaseSequence(this, MIN_VARIANT_LENGTH);
        }
        else
        {
            mRefBasesRepeatedTrimmed = formRefBaseSequence(MIN_VARIANT_LENGTH);
        }

        mRefBaseTrimLength = calcTrimmedRefBaseLength(this);
    }

    public String refBasesRepeatedTrimmed() { return mRefBasesRepeatedTrimmed; }
    public int refBaseTrimLength() { return mRefBaseTrimLength; }

    public List<RefSideSoftClip> refSideSoftClips() { return mRefSideSoftClips; }

    public boolean checkAddRefSideSoftClip(final Read read)
    {
        return RefSideSoftClip.checkAddRefSideSoftClip(mRefSideSoftClips, mJunction, read);
    }

    public List<RemoteRegion> remoteRegions() { return mRemoteRegions; }
    public void addRemoteRegions(final List<RemoteRegion> regions) { mRemoteRegions.addAll(regions); }

    public PhaseGroup phaseGroup() { return mPhaseGroup; }
    public void setPhaseGroup(final PhaseGroup phaseGroup) { mPhaseGroup = phaseGroup; }

    public PhaseSet phaseSet()
    {
        // can only be a part of one and could add a reference but for now is retrieved
        return mPhaseGroup != null ? mPhaseGroup.findPhaseSet(this) : null;
    }

    public AssemblyOutcome outcome() { return mOutcome; }

    public void setOutcome(final AssemblyOutcome outcome)
    {
        if(mOutcome != REMOTE_REGION && mOutcome != LOCAL_INDEL) // persist classification for now
            mOutcome = outcome;
    }

    public AlignmentOutcome alignmentOutcome() { return mAlignmentOutcome; }
    public void setAlignmentOutcome(final AlignmentOutcome outcome) { mAlignmentOutcome = outcome; }

    public void setAssemblyAlignmentInfo(final String info) { mAssemblyAlignmentInfo = info; }
    public String assemblyAlignmentInfo() { return mAssemblyAlignmentInfo != null ? mAssemblyAlignmentInfo : mJunction.coords(); }

    public void setReadIndices()
    {
        for(SupportRead read : mSupport)
        {
            if(read.type().isSplitSupport() || read.type() == EXTENSION)
            {
                // say assembly junc index = 100, junction index in read = 70, then read's start index in assembly is 100 - 70 = 30
                // for an extension read (where the junction position isn't in the read)
                int assemblyIndex = mJunctionIndex - read.junctionReadIndex();
                read.setJunctionAssemblyIndex(assemblyIndex);
            }
            else
            {
                // TODO: use read junction offset instead of calculating from scratch

                // set based on the relative positions
                if(isForwardJunction())
                {
                    int assemblyIndex = read.unclippedStart() - mRefBasePosition;
                    read.setJunctionAssemblyIndex(assemblyIndex);
                }
                else
                {
                    int junctionOffset = read.unclippedStart() - mJunction.Position;
                    int assemblyIndex = mJunctionIndex + junctionOffset;
                    read.setJunctionAssemblyIndex(assemblyIndex);
                }
            }
        }
    }

    public JunctionAssembly(
            final JunctionAssembly initialAssembly, final RefSideSoftClip refSideSoftClip, int refBaseLength,
            final List<SupportRead> initialSupport)
    {
        // build a junction assembly from an initial junction where the ref bases are truncated due to branching (likely short TI)
        mJunction = initialAssembly.junction();

        // copy the initial assembly's extension bases and ref bases up to the ref-side soft clip
        int extensionLength = initialAssembly.extensionLength();
        int newBaseLength = extensionLength + refBaseLength;
        int initialBaseLength = initialAssembly.baseLength();
        int baseLengthDiff = initialBaseLength - newBaseLength;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int assemblyIndexOffset = 0;

        if(mJunction.isForward())
        {
            assemblyIndexOffset = baseLengthDiff;
            mJunctionIndex = newBaseLength - extensionLength - 1;
        }
        else
        {
            mJunctionIndex = initialAssembly.junctionIndex();
        }

        mRefBasePosition = refSideSoftClip.Position;

        int assemblyIndex = assemblyIndexOffset;
        for(int i = 0; i < mBases.length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
            {
                mBases[i] = 0;
                mBaseQuals[i] = 0;
            }
            else
            {
                if(assemblyIndex >= initialAssembly.bases().length)
                    break;

                mBases[i] = initialAssembly.bases()[assemblyIndex];
                mBaseQuals[i] = initialAssembly.baseQuals()[assemblyIndex];
            }
        }

        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();

        mRepeatInfo = Lists.newArrayList();
        mRefBaseIndels = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = "";
        mRefBaseTrimLength = 0;
        mRefSideSoftClips = Lists.newArrayList(refSideSoftClip);
        mRemoteRegions = Lists.newArrayList();
        mMergedAssemblies = 0;
        mMismatchReadCount = 0;
        mPhaseGroup = null;

        mStats = new AssemblyStats();

        SupportRead initialRead = null;

        for(SupportRead support : initialSupport)
        {
            mSupport.add(support);

            mStats.addRead(support, mJunction, null);

            if(support.id().equals(initialAssembly.initialReadId()))
                initialRead = support;
        }

        mInitialReadId = initialRead != null ? initialRead.id() : (!mSupport.isEmpty() ? mSupport.get(0).id() : "");
        mIndelCoords = initialAssembly.indelCoords();
        mOutcome = UNSET;
        mAlignmentOutcome = NO_SET;
    }

    public void addCandidateSupport(final Read read)
    {
        mCandidateSupport.add(read);
        ++mStats.CandidateSupportCount;
    }

    public List<Read> candidateSupport() { return mCandidateSupport; }
    public void clearCandidateSupport() { mCandidateSupport.clear(); }

    public void addUnmappedRead(final Read read) { mUnmappedCandidates.add(read); }
    public List<Read> unmappedReads() { return mUnmappedCandidates; }

    public String toString()
    {
        return format("junc(%s) coords(extLen=%d refLen=%d refBasePos=%d len=%d juncIndex=%d) support(%d) mismatches(%d)",
                mJunction.coords(), extensionLength(), refBaseLength(), refBasePosition(), baseLength(), mJunctionIndex,
                mSupport.size(), mMismatchReadCount);
    }

    public String formFullSequence() { return formJunctionSequence(refBaseLength()); }

    public String formJunctionSequence()
    {
        return formJunctionSequence(0);
    }

    public String formJunctionSequence(final int includeRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mJunction.isForward())
        {
            seqIndexStart = max(0, mJunctionIndex + 1 - includeRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionIndex - 1 + includeRefBaseCount, mBases.length - 1);
        }

        return formSequence(seqIndexStart, seqIndexEnd);
    }

    public String formRefBaseSequence() { return formRefBaseSequence(-1); }

    public String formRefBaseSequence(int maxRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mJunction.isForward())
        {
            seqIndexEnd = mJunctionIndex;

            if(maxRefBaseCount > 0)
                seqIndexStart = seqIndexEnd - maxRefBaseCount + 1;
            else
                seqIndexStart = 0;
        }
        else
        {
            seqIndexStart = mJunctionIndex;

            if(maxRefBaseCount > 0)
                seqIndexEnd = seqIndexStart + maxRefBaseCount - 1;
            else
                seqIndexEnd = mBases.length - 1;
        }

        return formSequence(seqIndexStart, seqIndexEnd);
    }

    public String formSequence(int seqIndexStart, int seqIndexEnd)
    {
        StringBuilder sb = new StringBuilder();

        for(int index = max(seqIndexStart, 0); index <= min(seqIndexEnd, mBases.length - 1); ++index)
        {
            sb.append((char)mBases[index]);
        }

        return sb.toString();
    }

    // convenience
    public int minAlignedPosition() { return mJunction.isForward() ? refBasePosition() : mJunction.Position; }
    public int maxAlignedPosition() { return mJunction.isReverse() ? refBasePosition() : mJunction.Position; }

    @VisibleForTesting
    public JunctionAssembly(final Junction junction, final byte[] bases, final byte[] quals, final int junctionIndex)
    {
        mJunction = junction;
        mInitialReadId = null;

        mJunctionIndex = junctionIndex;

        // pos = 20, index = 10, length = 21, min align = 10, max align = 30

        if(mJunction.isForward())
        {
            mRefBasePosition = mJunction.Position - junctionIndex;
        }
        else
        {
            mRefBasePosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }

        mBases = copyArray(bases);
        mBaseQuals = copyArray(quals);
        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = "";
        mRefBaseTrimLength = 0;
        mRemoteRegions = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mRefBaseIndels = Lists.newArrayList();
        mMergedAssemblies = 0;
        mOutcome = UNSET;
        mAlignmentOutcome = NO_SET;
        mMismatchReadCount = 0;
        mStats = new AssemblyStats();
        mIndelCoords = null;
    }

    @VisibleForTesting
    public void addJunctionRead(final Read read)
    {
        int junctionReadIndex = mJunction.Position - read.unclippedStart();

        SupportRead support = new SupportRead(read, JUNCTION, junctionReadIndex, read.basesLength(), 0);
        mSupport.add(support);
        mStats.addRead(support, mJunction, read);
    }
}
