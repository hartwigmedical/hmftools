package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Character.toLowerCase;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.compareSequences;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class UnmappedBaseExtender
{
    private final JunctionAssembly mJunctionAssembly;
    private final boolean mIsForwardJunction;

    private byte[] mBases;
    private byte[] mBaseQuals;
    private final List<RepeatInfo> mRepeats;
    private final List<SupportRead> mSupportReads;

    private String mExtensionBases;
    private final int mRefBaseCount;

    public UnmappedBaseExtender(final JunctionAssembly assembly)
    {
        mJunctionAssembly = assembly;

        mIsForwardJunction = assembly.junction().isForward();

        int extensionIndexStart, extensionIndexEnd;

        // include a subset of the junction ref bases from the assembly plus all extension bases
        mRefBaseCount = min(ASSEMBLY_LINK_OVERLAP_BASES, assembly.refBaseLength());

        if(mIsForwardJunction)
        {
            extensionIndexStart = assembly.junctionIndex() - mRefBaseCount + 1;
            extensionIndexEnd = assembly.baseLength() - 1;
        }
        else
        {
            extensionIndexStart = 0;
            extensionIndexEnd = assembly.junctionIndex() + mRefBaseCount - 1;
        }

        mBases = subsetArray(assembly.bases(), extensionIndexStart, extensionIndexEnd);
        mBaseQuals = subsetArray(assembly.baseQuals(), extensionIndexStart, extensionIndexEnd);
        mSupportReads = Lists.newArrayList();

        mExtensionBases = null;
        mRepeats = Lists.newArrayList();
    }

    public byte[] extensionBases()
    {
        if(mIsForwardJunction)
        {
            // 0-9 ref, 10-19 ext, ref count = 10, total = 20
            return subsetArray(mBases, mRefBaseCount, mBases.length - 1);
        }
        else
        {
            // 0-9 ext, 10-19 ref, ref count = 10, total = 20
            return subsetArray(mBases, 0, mBases.length - mRefBaseCount - 1);
        }
    }

    public int extensionBaseLength() { return mBases.length - mRefBaseCount; }

    public byte[] baseQualities()
    {
        if(mIsForwardJunction)
            return subsetArray(mBaseQuals, mRefBaseCount, mBaseQuals.length - 1);
        else
            return subsetArray(mBaseQuals, 0, mBaseQuals.length - mRefBaseCount - 1);
    }

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
                int readExtensionLength = readSequenceMatch.maxReadExtension();

                if(readExtensionLength > 0)
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

        // build out extension bases from these overlapping reads
        int newExtensionLength = readSequenceMatches.stream().mapToInt(x -> x.maxReadExtension()).max().orElse(0);

        int baseOffset = 0;

        if(newExtensionLength > mBases.length)
        {
            if(!mIsForwardJunction)
                baseOffset =newExtensionLength - mBases.length;

            extendBases(newExtensionLength);
        }

        // now add reads to extend the new bases
        int addedReads = 0;

        Collections.sort(readSequenceMatches);

        for(ReadSequenceMatch readSequenceMatch : readSequenceMatches)
        {
            if(addRead(readSequenceMatch, baseOffset))
            {
                ++addedReads;
            }
        }

        trimUnsetBases();

        if(addedReads > 0 && AssemblyConfig.AssemblyBuildDebug)
        {
            SV_LOGGER.debug("junc({}) added {} unmapped reads, new sequence {}",
                    mJunctionAssembly.junction().coords(), addedReads, new String(mBases));
        }

        return addedReads > 0;
    }

    private void extendBases(int newExtensionLength)
    {
        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int baseOffset = newExtensionLength - mBases.length;

        mBases = new byte[newExtensionLength];
        mBaseQuals = new byte[newExtensionLength];

        for(int i = 0; i < newExtensionLength; ++i)
        {
            if(mIsForwardJunction)
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

    private void trimUnsetBases()
    {
        boolean isForward = mJunctionAssembly.isForwardJunction();
        int currentIndex = isForward ? 0 : mBases.length - 1;
        int validLength = 0;

        while(currentIndex >= 0 && currentIndex < mBases.length)
        {
            if(mBases[currentIndex] == 0)
                break;

            ++validLength;
            currentIndex += isForward ? 1 : -1;
        }

        if(validLength == mBases.length)
            return;

        int reduction = mBases.length - validLength;
        int startIndex = isForward ? 0 : reduction;
        int endIndex = startIndex + validLength - 1;
        mBases = subsetArray(mBases, startIndex, endIndex);
        mBaseQuals = subsetArray(mBaseQuals, startIndex, endIndex);
    }

    private boolean addRead(final ReadSequenceMatch readSequenceMatch, int baseOffset)
    {
        int mismatchCount = 0;

        // add in new bases to the extension sequence
        int readIndexStart, readIndexEnd, junctionReadStartDistance, extBaseIndexStart;

        Read read = readSequenceMatch.Read;

        if(mIsForwardJunction)
        {
            readIndexStart = readSequenceMatch.ReadSeqStart;
            readIndexEnd = read.basesLength() - 1;
            extBaseIndexStart = readSequenceMatch.ExtensionBaseSeqStart;
            junctionReadStartDistance = -readSequenceMatch.ExtensionBaseSeqStart - readSequenceMatch.ReadSeqStart;

            junctionReadStartDistance += mRefBaseCount;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = readSequenceMatch.Overlap - 1;
            extBaseIndexStart = baseOffset + readSequenceMatch.ExtensionBaseSeqStart - readSequenceMatch.ReadSeqStart;

            if(extBaseIndexStart < 0)
            {
                readIndexStart += abs(extBaseIndexStart);
                extBaseIndexStart = 0;
            }

            // say 0-9 is extension and junction index = 10
            // if extBaseIndex (read start vs extension bases) is 4, then junc index in read is
            junctionReadStartDistance = mBases.length - readSequenceMatch.ExtensionBaseSeqStart - readSequenceMatch.ReadSeqStart;

            // add the extension length since the read sequence match was prior to extension
            junctionReadStartDistance += baseOffset;

            junctionReadStartDistance -= mRefBaseCount;
        }

        boolean reverseBases = reverseReadBases(read);
        int readBaseLength = read.basesLength();

        int extBaseIndex = extBaseIndexStart;

        // first pass checks for mismatches only
        for(int i = readIndexStart; i <= readIndexEnd; ++i, ++extBaseIndex)
        {
            int readIndex = reverseBases ? readBaseLength - 1 - i : i;
            byte base = read.getBases()[readIndex];

            if(reverseBases)
                base = swapDnaBase(base);

            byte qual = read.getBaseQuality()[readIndex];

            if(extBaseIndex >= mBases.length)
                break;

            if(mBases[extBaseIndex] == 0)
                continue;

            if(mBases[extBaseIndex] == base || belowMinQual(qual) || belowMinQual(mBaseQuals[extBaseIndex]))
                continue;

            ++mismatchCount;
        }

        int readBaseOverlap = readIndexEnd - readIndexStart + 1;
        int permittedMismatches = permittedReadMismatches(readBaseOverlap);

        if(mismatchCount > permittedMismatches)
            return false;

        mismatchCount = 0;
        extBaseIndex = extBaseIndexStart;

        for(int i = readIndexStart; i <= readIndexEnd; ++i, ++extBaseIndex)
        {
            int readIndex = reverseBases ? readBaseLength - 1 - i : i;
            byte base = read.getBases()[readIndex];

            if(reverseBases)
                base = swapDnaBase(base);

            byte qual = read.getBaseQuality()[readIndex];

            if(extBaseIndex >= mBases.length)
                break;

            if(mBases[extBaseIndex] == 0)
            {
                mBases[extBaseIndex] = base;
                mBaseQuals[extBaseIndex] = qual;
            }
            else
            {
                if(mBases[extBaseIndex] == base || belowMinQual(qual))
                {
                    if((int)qual > (int)mBaseQuals[extBaseIndex])
                        mBaseQuals[extBaseIndex] = qual;
                }
                else if(belowMinQual(mBaseQuals[extBaseIndex]))
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

        if(AssemblyConfig.AssemblyBuildDebug)
        {
            char[] overlapReadBases = new char[readBaseOverlap];
            int baseIndex = 0;

            for(int i = readIndexStart; i <= readIndexEnd; ++i, ++extBaseIndex)
            {
                int readIndex = reverseBases ? readBaseLength - 1 - i : i;
                byte base = read.getBases()[readIndex];

                if(reverseBases)
                    base = swapDnaBase(base);

                byte qual = read.getBaseQuality()[readIndex];
                char baseChr = (char)base;

                overlapReadBases[baseIndex++] = belowMinQual(qual) ? toLowerCase(baseChr) : baseChr;
            }

            SV_LOGGER.debug("junc({}) added read({}) overlap({} index={}-{}), readBases: {}",
                    mJunctionAssembly.junction().coords(), readSequenceMatch.Read.id(), readBaseOverlap, readIndexStart, readIndexEnd,
                    new String(overlapReadBases));

            SV_LOGGER.debug("junc({}) added read({}) new sequence {}",
                    mJunctionAssembly.junction().coords(), readSequenceMatch.Read.id(), new String(mBases));
        }

        int matchedCount = readBaseOverlap - mismatchCount;
        SupportRead supportRead = new SupportRead(read, SupportType.EXTENSION, junctionReadStartDistance, matchedCount, mismatchCount);
        mSupportReads.add(supportRead);
        return true;
    }

    private class ReadSequenceMatch implements Comparable<ReadSequenceMatch>
    {
        public final Read Read;
        public final int ReadSeqStart;
        public final int ExtensionBaseSeqStart;
        public final int Overlap;
        public final int Mismatches;

        public ReadSequenceMatch(final Read read, final int readSeqStart, final int extensionBaseSeqStart, final int overlap, int mismatches)
        {
            Read = read;
            ReadSeqStart = readSeqStart;
            ExtensionBaseSeqStart = extensionBaseSeqStart;
            Overlap = overlap;
            Mismatches = mismatches;
        }

        public int maxReadExtension()
        {
            // say extension bases = 0-99 (100 in length), ext base match start = 50, corresponds to read index of 0, and read length = 151
            // then max read extension = 50 + 151 - 0 - 1 = 200
            if(mIsForwardJunction)
                return ExtensionBaseSeqStart + Read.basesLength() - ReadSeqStart;
            else
                return mBases.length - ExtensionBaseSeqStart + ReadSeqStart;
        }

        @Override
        public int compareTo(final ReadSequenceMatch other)
        {
            // longer overlap followed by lower mismatches
            if(Overlap != other.Overlap)
                return Overlap > other.Overlap ? -1 : 1;

            if(Mismatches != other.Mismatches)
                return Mismatches < other.Mismatches ? -1 : 1;

            return 0;
        }

        public String toString()
        {
            return format("id(%s) index(read=%d ext=%d) overlap(%d) mismatches(%d)",
                    Read.id(), ReadSeqStart, ExtensionBaseSeqStart, Overlap, Mismatches);
        }
    }

    private ReadSequenceMatch findReadSequenceOverlap(final Read read)
    {
        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;

        // apply qual trimming before a read is used
        ReadAdjustments.trimLowQualBases(read);

        int readBaseCount = read.getBases().length;

        boolean reverseBases = reverseReadBases(read);
        byte[] readBases = reverseBases ? reverseComplementBases(read.getBases()) : read.getBases();

        byte[] readBaseQuals = null;

        for(int readIndex = 0; readIndex + subsequenceLength < read.basesLength(); readIndex += subsequenceLength)
        {
            String readSeqBases = new String(readBases, readIndex, subsequenceLength);
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

            int permittedMismatches = permittedReadMismatches(totalOverlap);

            if(readBaseQuals == null)
                readBaseQuals = reverseBases ? reverseArray(read.getBaseQuality()) : read.getBaseQuality();

            int mismatchCount = compareSequences(
                    mBases, mBaseQuals, extBaseIndexStart, extBaseIndexEnd, mRepeats,
                    readBases, readBaseQuals, readIndexStart, readIndexEnd, Collections.emptyList(), permittedMismatches);

            if(mismatchCount > permittedMismatches)
                continue;

            return new ReadSequenceMatch(read, readIndexStart, extBaseIndexStart, totalOverlap, mismatchCount);
        }

        return null;
    }

    private int permittedReadMismatches(int readBaseOverlap)
    {
        // more lenient than normal read comparisons
        return mismatchesPerComparisonLength(readBaseOverlap) + 1;
    }

    private boolean reverseReadBases(final Read read) { return read.orientation() == read.mateOrientation(); }

    public String toString()
    {
        return format("junc(%s) extLen(%d) support(%d)",
                mJunctionAssembly.junction().coords(), mBases.length - mRefBaseCount, mSupportReads.size());
    }
}
