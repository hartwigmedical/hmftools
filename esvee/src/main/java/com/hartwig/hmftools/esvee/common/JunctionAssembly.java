package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.findUnsetBases;
import static com.hartwig.hmftools.esvee.common.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.common.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.read.ReadUtils.copyArray;

import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.filters.FilterType;
import com.hartwig.hmftools.esvee.read.Read;

public class JunctionAssembly
{
    private int mAssemblyId;
    private final Junction mJunction;
    private final Read mInitialRead;

    private int mJunctionSequenceIndex; // position of the junction in the read bases

    // aligned position on the extension side is an inferred position only
    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<AssemblySupport> mSupport;
    private final List<AssemblySupport> mCandidateSupport;
    private final List<RefSideSoftClip> mRefSideSoftClips;

    private final SequenceMismatches mSequenceMismatches;
    private final List<RepeatInfo> mRepeatInfo;

    private final List<RemoteRegion> mRemoteRegions;

    private PhaseGroup mPhaseGroup;
    private AssemblyOutcome mOutcome;

    private final List<JunctionAssembly> mBranchedAssemblies;

    // info only
    private int mMergedAssemblies;

    private final Set<FilterType> mFilters;

    public JunctionAssembly(
            final Junction junction, final Read read, final int maxExtensionDistance,
            final int minAlignedPosition, final int maxAlignedPosition)
    {
        mJunction = junction;
        mInitialRead = read;

        mJunctionSequenceIndex = junction.isForward() ? 0 : maxExtensionDistance;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        int initialAssemblyLength = maxExtensionDistance + 1;
        mBases = new byte[initialAssemblyLength];

        mBaseQuals = new byte[initialAssemblyLength];
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mBranchedAssemblies = Lists.newArrayList();
        mMergedAssemblies = 0;
        mPhaseGroup = null;
        mOutcome = UNSET;
        mFilters = Sets.newHashSet();

        addInitialRead(read);
    }

    public void setId(int id) { mAssemblyId = id; }
    public int id() { return mAssemblyId; }

    public Read initialRead() { return mInitialRead; }
    public Junction junction() { return mJunction; }
    public boolean isForwardJunction() { return mJunction.isForward(); }

    public boolean indel() { return mJunction.IndelBased; }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionSequenceIndex; };

    // eg 21 bases, junction index at 10 (so 0-9 = 10 before, 11-20 = 10 after), note: doesn't count the junction base
    public int lowerDistanceFromJunction() { return mJunctionSequenceIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionSequenceIndex - 1; };

    public int refBaseLength() { return (mJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction()) + 1; }

    // base count beyond the junction
    public int extensionLength() { return mJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    // purely for informational purposes at this stage - since past the junction is just the soft-clip length, and on the
    // reference side doesn't take into account any INDELs
    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int refBasePosition() { return mJunction.isForward() ? mMinAlignedPosition : mMaxAlignedPosition; }
    public int baseLength() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public SequenceMismatches mismatches() { return mSequenceMismatches; }

    private void addInitialRead(final Read read)
    {
        for(int i = 0; i < mBases.length; ++i)
        {
            mBases[i] = 0;
            mBaseQuals[i] = 0;
        }

        addReadSupport(read, mJunction.IndelBased ? INDEL : JUNCTION, false);
    }

    public boolean checkReadMatches(final Read read, int permittedMismatches)
    {
        int mismatchCount = 0;

        int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);
        int[] readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, false);
        int assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], false);

        if(assemblyIndex < 0)
        {
            // SV_LOGGER.debug("readCoords({}) invalid for assembly({}) read({})", readCoords, toString(), read.toString());
            return false;
        }

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length) // CHECK: similar to issue with INDELs in addRead() ??
                break;

            if(!basesMatch(
                    read.getBases()[i], mBases[assemblyIndex], read.getBaseQuality()[i], mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return true;
   }

    public void addReadSupport(final Read read, final SupportType type, boolean registerMismatches)
    {
        addRead(read, type, registerMismatches, null);
    }

    public void addJunctionRead(final Read read, boolean registerMismatches) { addReadSupport(read, JUNCTION, registerMismatches); }

    public void extendJunctionReadSupport(final Read read, final AssemblySupport existingSupport)
    {
        addRead(read, existingSupport.type(), true, existingSupport);
    }

    private void addRead(final Read read, final SupportType type, boolean registerMismatches, @Nullable final AssemblySupport existingSupport)
    {
        int mismatchCount = 0;
        int readJunctionIndex;
        int[] readIndexRange;
        int assemblyIndex;

        boolean byReference = existingSupport != null;

        if(type == JUNCTION || type == INDEL)
        {
            if(existingSupport == null)
            {
                readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);
            }
            else
            {
                readJunctionIndex = existingSupport.junctionReadIndex();
            }

            readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, byReference);

            assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], byReference);
        }
        else
        {
            readIndexRange = new int[] { 0, read.getBases().length - 1 }; // take the whole read for ref-side reads
            readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, false);

            if(mJunction.isForward())
            {
                assemblyIndex = read.unclippedStart() - mMinAlignedPosition;
            }
            else
            {
                assemblyIndex = mJunctionSequenceIndex + read.unclippedStart() - mJunction.Position;
            }
        }

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length)
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
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                }
                else
                {
                    ++mismatchCount;

                    if(registerMismatches)
                        mSequenceMismatches.add(assemblyIndex, base, read, qual);
                }
            }
        }

        if(existingSupport == null)
        {
            mSupport.add(new AssemblySupport(read, type, assemblyIndex, readJunctionIndex, readIndexRange, mismatchCount));
        }
        else
        {
            int[] existingReadRange = existingSupport.readIndexRange();

            existingSupport.setReadIndexRange(
                    min(existingReadRange[SE_START], readIndexRange[SE_START]),
                    max(existingReadRange[SE_END], readIndexRange[SE_END]));

            existingSupport.setReferenceMismatches(mismatchCount);
        }
    }

    public void extendBases(
            int refBaseExtensionDistance, int newMinAlignedPosition, int newMaxAlignedPosition, final RefBaseAssembly refBaseAssembly)
    {
        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int newBaseLength = mBases.length + refBaseExtensionDistance;

        boolean isForwardJunction = mJunction.isForward();

        if(isForwardJunction)
            mJunctionSequenceIndex += refBaseExtensionDistance;

        int baseOffset = isForwardJunction ? refBaseExtensionDistance : 0;

        int refBaseOffset = refBaseAssembly != null && isForwardJunction ? 0 : extensionLength();

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        mMinAlignedPosition = newMinAlignedPosition;
        mMaxAlignedPosition = newMaxAlignedPosition;

        // forward junction: 0 up to base offset - 1 -> set to zero, base offset to new length -> copy bases
        // reverse junction: 0 up to base offset - 1 -> copy bases, base offset to new length -> set to zero

        for(int i = 0; i < newBaseLength; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    if(refBaseAssembly != null && i < refBaseAssembly.bases().length)
                        mBases[i] = refBaseAssembly.bases()[i];
                    else
                        mBases[i] = 0;

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
                    int refBaseIndex = i - refBaseOffset;
                    if(refBaseAssembly != null && refBaseIndex < refBaseAssembly.bases().length)
                        mBases[i] = refBaseAssembly.bases()[refBaseIndex];
                    else
                        mBases[i] = 0;

                    mBaseQuals[i] = 0;
                }
            }
        }
    }

    public void mergeRefBaseAssembly(final RefBaseAssembly refBaseAssembly)
    {
        // find the longest length of aligned reference bases extending back from the junction
        int newRefBaseCount = refBaseAssembly.validRefBaseLength();
        int existingRefBaseCount = refBaseLength();

        if(newRefBaseCount > existingRefBaseCount)
        {
            if(isForwardJunction())
                mMinAlignedPosition = refBaseAssembly.minAlignedPosition();
            else
                mMaxAlignedPosition = refBaseAssembly.maxAlignedPosition();

            int refBaseExtension = newRefBaseCount - existingRefBaseCount;

            extendBases(refBaseExtension, mMinAlignedPosition, mMaxAlignedPosition, refBaseAssembly);
        }

        for(AssemblySupport support : refBaseAssembly.support())
        {
            addRead(support.read(), support.type(), true, null);
        }

        // use the ref assembly to fill in any missing bases

        List<int[]> emptyBaseRanges = findUnsetBases(mBases);

        if(!emptyBaseRanges.isEmpty())
        {
            SV_LOGGER.debug("assembly({}) refBases(existing={} new={}) empty ranges: {}",
                    refBaseAssembly, existingRefBaseCount, newRefBaseCount, emptyBaseRanges);
        }
    }

    private int[] getReadIndexCoordinates(final Read read, final int readJunctionIndex, boolean byReferenceBases)
    {
        // using the position of the junction within the read coordinates, gets the read index range either from the junction into
        // soft-clip / junction bases, or from the junction, but not including it, back into reference bases
        int[] readIndexCoords = {0, 0};

        if(byReferenceBases)
        {
            if(mJunction.isForward())
            {
                readIndexCoords[0] = 0;
                readIndexCoords[1] = readJunctionIndex - 1; // the base at the junction will have already been set
            }
            else
            {
                readIndexCoords[0] = readJunctionIndex + 1;
                readIndexCoords[1] = read.basesLength() - 1;
            }
        }
        else
        {
            if(mJunction.isForward())
            {
                readIndexCoords[0] = readJunctionIndex;
                readIndexCoords[1] = read.basesLength() - 1;
            }
            else
            {
                readIndexCoords[0] = 0;
                readIndexCoords[1] = readJunctionIndex;
            }
        }

        return readIndexCoords;
    }

    private int getReadAssemblyStartIndex(final int readJunctionIndex, final int readStartIndex, final boolean byReferenceBases)
    {
        if(byReferenceBases)
        {
            if(mJunction.isForward())
                return mJunctionSequenceIndex - (readJunctionIndex - readStartIndex);
            else
                return mJunctionSequenceIndex + 1;
        }
        else
        {
            if(mJunction.isForward())
                return mJunctionSequenceIndex; // ie bases starting from the junction
            else
                return mJunctionSequenceIndex - (readJunctionIndex - readStartIndex);
        }
    }

    public void removeSupportRead(final Read read)
    {
        for(int i = 0; i < mSupport.size(); ++i)
        {
            if(mSupport.get(i).read() == read)
            {
                mSupport.remove(i);
                return;
            }
        }
    }

    public void checkAddReadSupport(final JunctionAssembly other)
    {
        for(AssemblySupport support : other.support())
        {
            if(hasReadSupport(support.read()))
                continue;

            addJunctionRead(support.read(), false);
        }
    }

    public boolean hasReadSupport(final Read read)
    {
        return read != null && mSupport.stream().anyMatch(x -> x.read() == read);
    }

    // caching repeat info needs careful consideration since any extension of ref bases invalidates the values,
    // at least for +ve orientation assemblies
    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
            mRepeatInfo.addAll(repeats);
    }

    public List<RefSideSoftClip> refSideSoftClips() { return mRefSideSoftClips; }

    public boolean checkAddRefSideSoftClip(final Read read)
    {
        return RefSideSoftClip.checkAddRefSideSoftClip(mRefSideSoftClips, mJunction, read);
    }

    public List<RemoteRegion> remoteRegions() { return mRemoteRegions; }
    public void addRemoteRegions(final List<RemoteRegion> regions) { mRemoteRegions.addAll(regions); }

    public PhaseGroup phaseGroup() { return mPhaseGroup; }
    public void setPhaseGroup(final PhaseGroup phaseGroup) { mPhaseGroup = phaseGroup; }

    public AssemblyOutcome outcome() { return mOutcome; }
    public void setOutcome(final AssemblyOutcome outcome) { mOutcome = outcome; }

    public void addBranchedAssembly(final JunctionAssembly assembly)
    {
        mBranchedAssemblies.add(assembly);
    }

    public List<JunctionAssembly> branchedAssemblies() { return mBranchedAssemblies; }
    public boolean hasBranchedAssembly(final JunctionAssembly other) { return mBranchedAssemblies.contains(other); }

    public JunctionAssembly(
            final JunctionAssembly initialAssembly, final RefSideSoftClip refSideSoftClip,
            final List<AssemblySupport> initialSupport, final Set<Read> excludedReads)
    {
        // build a junction assembly from an initial junction where the ref bases are truncated due to branching (likely short TI)
        mJunction = initialAssembly.junction();

        // copy the initial assembly's extension bases and ref bases up to the ref-side soft clip
        int extensionLength = initialAssembly.extensionLength();
        int newBaseLength = extensionLength + abs(mJunction.Position - refSideSoftClip.Position) + 1;
        int initialBaseLength = initialAssembly.baseLength();

        int baseLengthDiff = initialBaseLength - newBaseLength;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int assemblyIndexOffset = 0;

        if(mJunction.isForward())
        {
            assemblyIndexOffset = baseLengthDiff;
            mJunctionSequenceIndex = newBaseLength - extensionLength - 1;
            mMinAlignedPosition = refSideSoftClip.Position;
            mMaxAlignedPosition = initialAssembly.maxAlignedPosition();
        }
        else
        {
            mJunctionSequenceIndex = initialAssembly.junctionIndex();
            mMinAlignedPosition = initialAssembly.minAlignedPosition();
            mMaxAlignedPosition = refSideSoftClip.Position;
        }

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

        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();

        mRepeatInfo = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList(refSideSoftClip);
        mRemoteRegions = Lists.newArrayList();
        mBranchedAssemblies = Lists.newArrayList();
        mFilters = Sets.newHashSet();
        mMergedAssemblies = 0;
        mPhaseGroup = null;

        Read initialRead = null;

        for(AssemblySupport support : initialSupport)
        {
            if(excludedReads.contains(support.read()))
                continue;

            // any reference mismatches aren't take since they were against the orginal assembly's ref bases
            mSupport.add(new AssemblySupport(
                    support.read(), support.type(), support.assemblyIndex(), support.junctionReadIndex(), support.readIndexRange(), support.junctionMismatches()));

            if(support.read() == initialAssembly.initialRead())
                initialRead = support.read();
        }

        if(initialRead == null && !mSupport.isEmpty())
            initialRead = mSupport.get(0).read();

        mInitialRead = initialRead;
        mOutcome = initialAssembly.outcome();
    }

    public void addCandidateSupport(final Read read, final SupportType type)
    {
        mCandidateSupport.add(new AssemblySupport(read, type, 0, 0, new int[] {0, 0}, 0));
    }

    public List<AssemblySupport> candidateSupport() { return mCandidateSupport; }
    public void clearCandidateSupport() { mCandidateSupport.clear(); }

    public Set<FilterType> filters() { return mFilters; }
    public boolean passing() { return mFilters.isEmpty(); }
    public void addFilter(final FilterType filterType) { mFilters.add(filterType); }

    public String toString()
    {
        return format("junc(%s) coords(extLen=%d refBasePos=%d len=%d) juncIndex(%d) support(%d) mismatches(pos=%d all=%d)",
                mJunction, extensionLength(), refBasePosition(), baseLength(), mJunctionSequenceIndex,
                mSupport.size(), mSequenceMismatches.positionCount(), mSequenceMismatches.distinctBaseCount());
    }

    public boolean hasUnsetBases() { return !findUnsetBases(mBases).isEmpty(); }

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
            seqIndexStart = max(0, mJunctionSequenceIndex + 1 - includeRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionSequenceIndex - 1 + includeRefBaseCount, mBases.length - 1);
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
            seqIndexEnd = mJunctionSequenceIndex;

            if(maxRefBaseCount > 0)
                seqIndexStart = seqIndexEnd - maxRefBaseCount + 1;
            else
                seqIndexStart = 0;
        }
        else
        {
            seqIndexStart = mJunctionSequenceIndex;

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


    @VisibleForTesting
    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] quals, final int junctionIndex)
    {
        mJunction = junction;
        mInitialRead = null;

        mJunctionSequenceIndex = junctionIndex;

        // pos = 20, index = 10, length = 21, min align = 10, max align = 30

        if(mJunction.isForward())
        {
            mMinAlignedPosition = mJunction.Position - junctionIndex;
            mMaxAlignedPosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }
        else
        {
            mMinAlignedPosition = mJunction.Position - junctionIndex;
            mMaxAlignedPosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }

        mBases = copyArray(bases);
        mBaseQuals = copyArray(quals);
        mSequenceMismatches = new SequenceMismatches();
        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mBranchedAssemblies = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mFilters = Sets.newHashSet();
        mMergedAssemblies = 0;
        mOutcome = UNSET;
    }
}
