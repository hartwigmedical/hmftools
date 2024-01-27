package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.read.ReadUtils.copyArray;

import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class JunctionAssembly
{
    private final Junction mInitialJunction;
    private final Read mInitialRead;

    private int mJunctionSequenceIndex;
    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];
    private int mBaseQualTotals[];

    private final List<AssemblySupport> mSupport;

    private final SequenceMismatches mSequenceMismatches;
    private final List<RepeatInfo> mRepeatInfo;

    // info only
    private int mMergedAssemblies;

    public JunctionAssembly(final Junction initialJunction, final Read read, final int maxExtensionDistance,
            final int minAlignedPosition, final int maxAlignedPosition)
    {
        mInitialJunction = initialJunction;
        mInitialRead = read;

        mJunctionSequenceIndex = initialJunction.isForward() ? 0 : maxExtensionDistance;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        int initialAssemblyLength = maxExtensionDistance + 1;
        mBases = new byte[initialAssemblyLength];

        mBaseQuals = new byte[initialAssemblyLength];
        mBaseQualTotals = new int[initialAssemblyLength];
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mMergedAssemblies = 0;

        addInitialRead(read);
    }

    public Read initialRead() { return mInitialRead; }
    public Junction initialJunction() { return mInitialJunction; }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionSequenceIndex; };

    // eg 21 bases, junction index at 10 (so 0-9 before, 11-20 after)
    public int lowerDistanceFromJunction() { return mJunctionSequenceIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionSequenceIndex - 1; };
    public int refBaseLength() { return mInitialJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction(); }
    public int extensionLength() { return mInitialJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    // purely for informational purposes at this stage - since past the junction is just the soft-clip length, and on the
    // reference side doesn't take into account any INDELs
    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int length() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }
    public int[] baseQualTotals() { return mBaseQualTotals; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public SequenceMismatches mismatches() { return mSequenceMismatches; }
    public boolean hasMismatches() { return mSequenceMismatches.hasMismatches(); }

    private void addInitialRead(final Read read)
    {
        for(int i = 0; i < mBases.length; ++i)
        {
            mBases[i] = 0;
            mBaseQuals[i] = 0;
            mBaseQualTotals[i] = 0;
        }

        addRead(read, false);
    }

    public boolean checkReadMatches(final Read read, int permittedMismatches)
    {
        int mismatchCount = 0;

        int readJunctionIndex = read.getReadIndexAtReferencePosition(mInitialJunction.Position, true);
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

    public void addRead(final Read read, boolean registerMismatches)
    {
        addRead(read, registerMismatches, null);
    }

    private void addRead(final Read read, boolean registerMismatches, @Nullable final AssemblySupport existingSupport)
    {
        int mismatchCount = 0;

        int readJunctionIndex;

        boolean byReference = existingSupport != null;

        if(existingSupport == null)
        {
            readJunctionIndex = read.getReadIndexAtReferencePosition(mInitialJunction.Position, true);
        }
        else
        {
            readJunctionIndex = existingSupport.junctionReadIndex();
        }

        int[] readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, byReference);

        int assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], byReference);

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length)
            {
                // can happen with INDELs in the CIGAR, need to adjust the read index ranges or iterate by elements
                // SV_LOGGER.debug("i({}) readCoords({}) assembly({}) read({})", i, readCoords, toString(), read.toString());
                break;
            }

            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;
                mBaseQualTotals[assemblyIndex] = qual;
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    mBaseQualTotals[assemblyIndex] += qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                    mBaseQualTotals[assemblyIndex] += qual;
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
            mSupport.add(new AssemblySupport(read, assemblyIndex, readJunctionIndex, readIndexRange, mismatchCount));
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

    public void expandReferenceBases()
    {
        // find the long length of reference bases extending back from the junction
        int minAlignedPosition = mMinAlignedPosition;
        int maxAlignedPosition = mMaxAlignedPosition;

        boolean isForwardJunction = mInitialJunction.isForward();
        int junctionPosition = mInitialJunction.Position;
        int maxDistanceFromJunction = 0;

        AssemblySupport minNmSupport = null;

        for(AssemblySupport support : mSupport)
        {
            int readJunctionIndex = support.read().getReadIndexAtReferencePosition(junctionPosition, true);

            // for positive orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            // for negative orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 4
            int extensionDistance = isForwardJunction ? readJunctionIndex : support.read().basesLength() - readJunctionIndex - 1;

            maxDistanceFromJunction = max(maxDistanceFromJunction, extensionDistance);

            if(isForwardJunction)
            {
                minAlignedPosition = min(minAlignedPosition, support.read().unclippedStart());
            }
            else
            {
                maxAlignedPosition = max(maxAlignedPosition, support.read().unclippedEnd());
            }

            if(minNmSupport == null || support.read().numberOfEvents() < minNmSupport.read().numberOfEvents())
            {
                minNmSupport = support;
            }
        }

        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);
        int[] existingBaseTotals = copyArray(mBaseQualTotals);

        int newBaseLength = mBases.length + maxDistanceFromJunction;

        if(isForwardJunction)
            mJunctionSequenceIndex += maxDistanceFromJunction;

        int baseOffset = isForwardJunction ? maxDistanceFromJunction : 0;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];
        mBaseQualTotals = new int[newBaseLength];

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        // forward junction: 0 up to base offset - 1 -> set to zero, base offset to new length -> copy bases
        // reverse junction: 0 up to base offset - 1 -> copy bases, base offset to new length -> set to zero

        for(int i = 0; i < newBaseLength; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    mBases[i] = 0;
                    mBaseQuals[i] = 0;
                    mBaseQualTotals[i] = 0;
                }
                else
                {
                    mBases[i] = existingBases[i - baseOffset];
                    mBaseQuals[i] = existingQuals[i - baseOffset];
                    mBaseQualTotals[i] = existingBaseTotals[i - baseOffset];
                }
            }
            else
            {
                if(i < existingBases.length)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                    mBaseQualTotals[i] = existingBaseTotals[i];
                }
                else
                {
                    mBases[i] = 0;
                    mBaseQuals[i] = 0;
                    mBaseQualTotals[i] = 0;
                }
            }
        }

        // order by NM to favour the ref where possible
        if(minNmSupport != null)
        {
            addRead(minNmSupport.read(), false, minNmSupport);
        }

        for(AssemblySupport support : mSupport)
        {
            if(support == minNmSupport)
                continue;

            addRead(support.read(), true, support);
        }
    }

    private int[] getReadIndexCoordinates(final Read read, final int readJunctionIndex, boolean byReferenceBases)
    {
        // using the position of the junction within the read coordinates, gets the read index range either from the junction into
        // soft-clip / junction bases, or from the junction, but not including it, back into reference bases
        int[] readIndexCoords = {0, 0};

        if(byReferenceBases)
        {
            if(mInitialJunction.isForward())
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
            if(mInitialJunction.isForward())
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
            if(mInitialJunction.isForward())
                return mJunctionSequenceIndex - (readJunctionIndex - readStartIndex);
            else
                return mJunctionSequenceIndex + 1;
        }
        else
        {
            if(mInitialJunction.isForward())
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
            if(mSupport.stream().anyMatch(x -> x.read() == support.read()))
                continue;

            addRead(support.read(), false);
        }
    }

    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
            mRepeatInfo.addAll(repeats);
    }

    public String toString()
    {
        return format("junc(%s) range(%d - %d len=%d) juncIndex(%d) support(%d) mismatches(pos=%d all=%d)",
                mInitialJunction, mMinAlignedPosition, mMaxAlignedPosition, length(), mJunctionSequenceIndex,
                mSupport.size(), mSequenceMismatches.positionCount(), mSequenceMismatches.distinctBaseCount());
    }

    public String formSequence(final int maxRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mInitialJunction.isForward())
        {
            seqIndexStart = max(0, mJunctionSequenceIndex - maxRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionSequenceIndex + maxRefBaseCount, mBases.length - 1);
        }

        StringBuilder sb = new StringBuilder();

        for(int index = seqIndexStart; index <= seqIndexEnd; ++index)
        {
            sb.append((char)mBases[index]);
        }

        return sb.toString();
    }
}
