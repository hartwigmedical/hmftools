package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.read.ReadUtils.copyArray;

import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class JunctionAssembly
{
    private final Junction mInitialJunction;
    private final Read mInitialRead;

    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];
    private int mBaseQualTotals[];

    private final List<AssemblySupport> mSupport;

    private final SequenceMismatches mSequenceMismatches;

    public JunctionAssembly(final Junction initialJunction, final Read read, final int minAlignedPosition, final int maxAlignedPosition)
    {
        mInitialJunction = initialJunction;
        mInitialRead = read;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        int initialAssemblyLength = maxAlignedPosition - minAlignedPosition + 1;
        mBases = new byte[initialAssemblyLength];

        mBaseQuals = new byte[initialAssemblyLength];
        mBaseQualTotals = new int[initialAssemblyLength];
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();

        addInitialRead(read);
    }

    public Read initialRead() { return mInitialRead; }
    public Junction initialJunction() { return mInitialJunction; }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int length() { return mMaxAlignedPosition - mMinAlignedPosition + 1; }

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

        int junctionReadIndex = read.getReadIndexAtReferencePosition(mInitialJunction.Position, true);
        int assemblyIndex = getReadAssemblyStartIndex(read, junctionReadIndex, false);

        if(assemblyIndex < 0)
        {
            // SV_LOGGER.debug("readCoords({}) invalid for assembly({}) read({})", readCoords, toString(), read.toString());
            return false;
        }

        int[] readIndexRange = getReadIndexCoordinates(read, junctionReadIndex, false);

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(i >= mBases.length) // CHECK: similar to issue with INDELs in addRead() ??
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

   private int[] getReadIndexCoordinates(final Read read, final int junctionReadIndex, boolean byReferenceBases)
   {
       int[] readIndexCoords = {0, 0};

       if(byReferenceBases)
       {
           if(mInitialJunction.isForward())
           {
               readIndexCoords[0] = 0;
               readIndexCoords[1] = junctionReadIndex - 1; // the base at the junction will have already been set
           }
           else
           {
               readIndexCoords[0] = junctionReadIndex + 1;
               readIndexCoords[1] = read.basesLength() - 1;
           }
       }
       else
       {
           if(mInitialJunction.isForward())
           {
               readIndexCoords[0] = junctionReadIndex;
               readIndexCoords[1] = read.basesLength() - 1;
           }
           else
           {
               readIndexCoords[0] = 0;
               readIndexCoords[1] = junctionReadIndex;
           }
       }

       return readIndexCoords;
   }

   private int getReadAssemblyStartIndex(final Read read, final int junctionReadIndex, boolean byReferenceBases)
   {
       if(byReferenceBases)
       {
           if(mInitialJunction.isForward())
               return read.unclippedStart() - mMinAlignedPosition;
           else
               return mInitialJunction.position() - mMinAlignedPosition + 1;
       }
       else
       {
           if(mInitialJunction.isForward())
               return 0;
           else
               return read.unclippedStart() - mMinAlignedPosition;
       }
   }

    public void addRead(final Read read, boolean registerMismatches)
    {
        addRead(read, registerMismatches, null);
    }

    private void addRead(final Read read, boolean registerMismatches, @Nullable final AssemblySupport existingSupport)
    {
        int mismatchCount = 0;

        int junctionReadIndex;

        boolean byReference = existingSupport != null;

        if(existingSupport == null)
        {
            junctionReadIndex = read.getReadIndexAtReferencePosition(mInitialJunction.Position, true);
        }
        else
        {
            junctionReadIndex = existingSupport.junctionReadIndex();
        }

        int assemblyIndex = getReadAssemblyStartIndex(read, junctionReadIndex, byReference);
        int[] readIndexRange = getReadIndexCoordinates(read, junctionReadIndex, byReference);

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
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
                else if(registerMismatches)
                {
                    // register a mismatch
                    ++mismatchCount;

                    if(registerMismatches)
                        mSequenceMismatches.add(assemblyIndex, base, read, qual);
                }
            }
        }

        if(existingSupport == null)
        {
            mSupport.add(new AssemblySupport(read, assemblyIndex, junctionReadIndex, readIndexRange, mismatchCount));
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
        int minAlignedPosition = mMinAlignedPosition;
        int maxAlignedPosition = mMaxAlignedPosition;

        boolean isForwardJunction = mInitialJunction.isForward();

        AssemblySupport minNmSupport = null;

        for(AssemblySupport support : mSupport)
        {
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

        // consider adjusting the assembly indices to be zero based, but actually keeping them based around zero for the junction
        // may make more sense

        int newBaseLength = maxAlignedPosition - minAlignedPosition + 1;
        int baseOffset = isForwardJunction ? mMinAlignedPosition - minAlignedPosition : 0;

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

            addRead(support.read(), false, support);
        }
    }

    public String toString()
    {
        return format("junc(%s) range(%d - %d len=%d) support(%d) mismatches(pos=%d all=%d)",
                mInitialJunction, mMinAlignedPosition, mMaxAlignedPosition, length(),
                mSupport.size(), mSequenceMismatches.positionCount(), mSequenceMismatches.distinctBaseCount());
    }

    public String formSequence(final int positionStart, final int positionEnd)
    {
        if(!positionsWithin(positionStart, positionEnd, mMinAlignedPosition, mMaxAlignedPosition))
            return "";

        StringBuilder sb = new StringBuilder();

        for(int pos = positionStart; pos <= positionEnd; ++pos)
        {
            sb.append((char)mBases[pos - mMinAlignedPosition]);
        }

        return sb.toString();
    }

    public boolean matches(final JunctionAssembly other, int minOverlapBases, int maxMismatchCount)
    {
        int overlapMinIndex = max(mMinAlignedPosition, other.minAlignedPosition());
        int overlapMaxIndex = min(mMaxAlignedPosition, other.maxAlignedPosition());

        int overlapLength = overlapMaxIndex - overlapMinIndex + 1;
        if(overlapLength < minOverlapBases)
            return false;

        int offset = overlapMinIndex - minAlignedPosition();
        int otherOffset = overlapMinIndex - other.minAlignedPosition();

        int mismatchCount = 0;

        for(int i = 0; i < overlapLength; ++i)
        {
            if(mBases[i + offset] != other.bases()[i + otherOffset])
            {
                ++mismatchCount;

                if(mismatchCount > maxMismatchCount)
                    return false;
            }
        }

        return true;
    }


}
