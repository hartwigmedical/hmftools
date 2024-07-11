package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.compareSequences;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class UnmappedBaseExtender
{
    private final JunctionAssembly mJunctionAssembly;
    private final Orientation mJunctionOrientation;

    private byte[] mBases;
    private byte[] mBaseQuals;
    private final List<RepeatInfo> mRepeats;
    private final List<SupportRead> mSupportReads;

    private String mExtensionBases;

    public UnmappedBaseExtender(final JunctionAssembly assembly)
    {
        mJunctionAssembly = assembly;

        mJunctionOrientation = assembly.junction().Orient;

        int extensionIndexStart, extensionIndexEnd;

        if(mJunctionOrientation.isForward())
        {
            extensionIndexStart = assembly.junctionIndex();
            extensionIndexEnd = assembly.baseLength() - 1;
        }
        else
        {
            extensionIndexStart = 0;
            extensionIndexEnd = assembly.junctionIndex();
        }

        mBases = subsetArray(assembly.bases(), extensionIndexStart, extensionIndexEnd);
        mBaseQuals = subsetArray(assembly.baseQuals(), extensionIndexStart, extensionIndexEnd);
        mSupportReads = Lists.newArrayList();

        mExtensionBases = null;
        mRepeats = Lists.newArrayList();
    }

    public byte[] extensionBases() { return mBases; }
    public int extensionBaseLength() { return mBases.length - 1; }
    public byte[] baseQualities() { return mBaseQuals; }
    public List<SupportRead> supportReads() { return mSupportReads; }

    public void processReads(final List<Read> reads)
    {
        while(true)
        {
            if(!tryAddReads(reads))
                break;
        }
    }

    public boolean tryAddReads(final List<Read> reads)
    {
        mExtensionBases = new String(mBases);

        mRepeats.clear();
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);
        if(repeats != null)
            mRepeats.addAll(repeats);

        List<ReadSequenceMatch> readSequenceMatches = Lists.newArrayList();

        int readIndex = 0;
        while(readIndex < reads.size())
        {
            Read read = reads.get(readIndex);

            ReadSequenceMatch readSequenceMatch = findReadSequenceOverlap(read);

            if(readSequenceMatch != null)
            {
                readSequenceMatches.add(readSequenceMatch);
                reads.remove(readIndex);
            }
            else
            {
                ++readIndex;
            }
        }

        if(readSequenceMatches.isEmpty())
            return false;

        Collections.sort(readSequenceMatches, Comparator.comparingInt(x -> -x.Overlap));

        // build out extension bases from these overlapping reads
        int newExtensionLength = readSequenceMatches.stream().mapToInt(x -> x.maxReadExtension()).max().orElse(0);
        int baseOffset = newExtensionLength - mBases.length;
        extendBases(newExtensionLength);

        // now add reads to extend the new bases
        readSequenceMatches.forEach(x -> addRead(x, baseOffset));

        return true;
    }

    private void extendBases(int newExtensionLength)
    {
        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int baseOffset = newExtensionLength - mBases.length;

        mBases = new byte[newExtensionLength];
        mBaseQuals = new byte[newExtensionLength];

        // adjust the extension base index for existing unmapped reads
        if(mJunctionOrientation.isReverse())
        {
            for(SupportRead read : mSupportReads)
            {
                read.setJunctionAssemblyIndex(read.junctionAssemblyIndex() + baseOffset);
            }
        }

        for(int i = 0; i < newExtensionLength; ++i)
        {
            if(mJunctionOrientation.isForward())
            {
                if(i < existingBases.length)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                }
                else
                {
                    mBases[i] = 0;
                    mBaseQuals[i] = 0;
                }
            }
            else
            {
                if(i < baseOffset)
                {
                    mBases[i] = 0;
                    mBaseQuals[i] = 0;
                }
                else
                {
                    mBases[i] = existingBases[i - baseOffset];
                    mBaseQuals[i] = existingQuals[i - baseOffset];
                }
            }
        }
    }

    private void addRead(final ReadSequenceMatch readSequenceMatch, int baseOffset)
    {
        int mismatchCount = 0;

        // add in new bases to the extension sequence
        int readIndexStart, readIndexEnd, junctionIndexOffset, extBaseIndex;

        Read read = readSequenceMatch.Read;

        if(mJunctionOrientation.isForward())
        {
            readIndexStart = readSequenceMatch.ReadSeqStart + readSequenceMatch.Overlap;
            readIndexEnd = read.basesLength() - 1;
            extBaseIndex = readSequenceMatch.ExtensionBaseSeqStart + readSequenceMatch.Overlap;
            junctionIndexOffset = -readSequenceMatch.ExtensionBaseSeqStart;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = readSequenceMatch.ReadSeqStart - 1;
            extBaseIndex = baseOffset + readSequenceMatch.ExtensionBaseSeqStart - readSequenceMatch.ReadSeqStart;

            if(extBaseIndex < 0)
            {
                readIndexStart += abs(extBaseIndex);
                extBaseIndex = 0;
            }

            // say 0-9 is extension and junction index = 10
            // if extBaseIndex (read start vs extension bases) is 4, then junc index in read is
            junctionIndexOffset = mBases.length - extBaseIndex;
        }

        for(int i = readIndexStart; i <= readIndexEnd; ++i, ++extBaseIndex)
        {
            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(extBaseIndex >= mBases.length)
                break;

            if(mBases[extBaseIndex] == 0)
            {
                mBases[extBaseIndex] = base;
                mBaseQuals[extBaseIndex] = qual;
            }
            else
            {
                if(mBases[extBaseIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[extBaseIndex])
                        mBaseQuals[extBaseIndex] = qual;
                }
                else if(mBaseQuals[extBaseIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[extBaseIndex] = base;
                    mBaseQuals[extBaseIndex] = qual;
                }
                else
                {
                    ++mismatchCount;
                }
            }
        }

        int readBaseOverlap = readIndexEnd - readIndexStart + 1;
        int permittedMismatches = mismatchesPerComparisonLength(readBaseOverlap);

        if(mismatchCount <= permittedMismatches)
        {
            int matchedCount = readBaseOverlap - mismatchCount;
            SupportRead supportRead = new SupportRead(read, SupportType.EXTENSION, junctionIndexOffset, matchedCount, mismatchCount);
            supportRead.setJunctionAssemblyIndex(extBaseIndex);
            mSupportReads.add(supportRead);
        }
    }

    private class ReadSequenceMatch
    {
        public final Read Read;
        public final int ReadSeqStart;
        public final int ExtensionBaseSeqStart;
        public final int Overlap;

        public ReadSequenceMatch(final Read read, final int readSeqStart, final int extensionBaseSeqStart, final int overlap)
        {
            Read = read;
            ReadSeqStart = readSeqStart;
            ExtensionBaseSeqStart = extensionBaseSeqStart;
            Overlap = overlap;
        }

        public int maxReadExtension()
        {
            // say extension bases = 0-99 (100 in length), ext base match start = 50, corresponds to read index of 0, and read length = 151
            // then max read extension = 50 + 151 - 0 - 1 = 200
            if(mJunctionOrientation.isForward())
                return ExtensionBaseSeqStart + Read.basesLength() - ReadSeqStart - 1;
            else
                return (mBases.length - ExtensionBaseSeqStart - 1) + ReadSeqStart;
        }

        public String toString()
        {
            return format("id(%s) index(read=%d ext=%d) overlap(%d)", Read.id(), ReadSeqStart, ExtensionBaseSeqStart, Overlap);
        }
    }

    private ReadSequenceMatch findReadSequenceOverlap(final Read read)
    {
        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;

        // String readBases = new String(read.getBases());
        int readBaseCount = read.getBases().length;

        for(int readIndex = 0; readIndex + subsequenceLength < read.basesLength(); readIndex += subsequenceLength)
        {
            String readSeqBases = new String(read.getBases(), readIndex, subsequenceLength);
            int extBaseMatchIndex = mExtensionBases.indexOf(readSeqBases);

            if(extBaseMatchIndex < 0)
                continue;

            // try to build out a longer sequence match from this point
            int readLowerBaseCount = readIndex;
            int readUpperBaseCount = readBaseCount - readIndex;
            int extBaseLowerBaseCount = extBaseMatchIndex;
            int extBaseUpperBaseCount = mBases.length - extBaseMatchIndex;

            int lowerOverlap = min(readLowerBaseCount, extBaseLowerBaseCount);
            int upperOverlap = min(readUpperBaseCount, extBaseUpperBaseCount);

            int totalOverlap = lowerOverlap + upperOverlap;

            if(totalOverlap < ASSEMBLY_LINK_OVERLAP_BASES)
                continue;

            int readIndexStart = readIndex - lowerOverlap;
            int readIndexEnd = readIndex + upperOverlap - 1;
            int extBaseIndexStart = extBaseMatchIndex - lowerOverlap;
            int extBaseIndexEnd = extBaseMatchIndex + upperOverlap - 1;

            int permittedMismatches = mismatchesPerComparisonLength(totalOverlap);

            int mismatchCount = compareSequences(
                    mBases, mBaseQuals, extBaseIndexStart, extBaseIndexEnd, mRepeats,
                    read.getBases(), read.getBaseQuality(), readIndexStart, readIndexEnd, Collections.emptyList(), permittedMismatches);

            if(mismatchCount > permittedMismatches)
                continue;

            return new ReadSequenceMatch(read, readIndexStart, extBaseIndexStart, totalOverlap);
        }

        return null;
    }
}
