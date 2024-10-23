package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_REF_READ_MIN_SOFT_CLIP;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class RefBaseSeqBuilder
{
    private final JunctionAssembly mAssembly;
    private final List<RefReadParseState> mReads;

    private final boolean mIsForward;
    private final int mJunctionPosition;

    private byte[] mBases; // starts with the junction ref base
    private byte[] mBaseQuals;
    private final List<CigarElement> mCigarElements;
    private int mRefBasePosition;

    public RefBaseSeqBuilder(final JunctionAssembly assembly)
    {
        mAssembly = assembly;
        mIsForward = mAssembly.junction().isForward();
        mJunctionPosition = mAssembly.junction().Position;
        mRefBasePosition = mJunctionPosition;

        int maxReadBaseLength = 0;

        // all reads will be loaded even if not used, to make updating ref bases more efficient back into the junction assembly
        mReads = Lists.newArrayListWithCapacity(assembly.supportCount());

        for(SupportRead support : assembly.support())
        {
            Read read = support.cachedRead();
            int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunctionPosition, true);

            boolean hasValidJunctionOverlap;

            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null && indelCoords.matchesJunction(mJunctionPosition, assembly.junction().Orient))
            {
                hasValidJunctionOverlap = true;
            }
            else
            {
                if(mIsForward)
                {
                    // junction reads must overlap the junction by 3+ bases to extend the ref sequence
                    hasValidJunctionOverlap = mAssembly.discordantOnly() || read.isRightClipped();
                    hasValidJunctionOverlap &= read.unclippedEnd() - mJunctionPosition >= ASSEMBLY_REF_READ_MIN_SOFT_CLIP;
                }
                else
                {
                    hasValidJunctionOverlap = mAssembly.discordantOnly() || read.isLeftClipped();
                    hasValidJunctionOverlap &= mJunctionPosition - read.unclippedStart() >= ASSEMBLY_REF_READ_MIN_SOFT_CLIP;
                }
            }

            int readRefBaseLength = readRefBaseLength(read, readJunctionIndex, mIsForward);

            RefReadParseState readState = new RefReadParseState(mIsForward, read, readJunctionIndex, readRefBaseLength);

            if(hasValidJunctionOverlap)
                maxReadBaseLength = max(maxReadBaseLength, readRefBaseLength);
            else
                readState.markInvalid();

            mReads.add(readState);
        }

        int baseLength = maxReadBaseLength + 1; // since the junction base itself is included (which is a ref base)

        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];
        mCigarElements = Lists.newArrayList();

        buildSequence();

        trimFinalSequence();
    }

    public byte[] bases() { return mBases; }
    public byte[] baseQualities() { return mBaseQuals; }
    public int refBaseLength() { return mBases.length; }
    public int refBasePosition() { return mRefBasePosition; }
    public List<CigarElement> cigarElements() { return mCigarElements; }
    public String cigarStr() { return CigarUtils.cigarElementsToStr(mCigarElements); }

    public static int readRefBaseLength(final Read read, int readJunctionIndex, boolean isForwardJunction)
    {
        // work out the aligned length including deletes and adding inserts
        int insertLength = read.cigarElements().stream().filter(x -> x.getOperator() == I).mapToInt(x -> x.getLength()).sum();

        // exclude soft clipped bases on the ref side
        int alignedLength = isForwardJunction ?
                max(readJunctionIndex - read.leftClipLength(), 0)
                : max(read.basesLength() - readJunctionIndex - 1 - read.rightClipLength(), 0);

        return alignedLength + insertLength;
    }

    public List<RefReadParseState> reads() { return mReads; }

    private static final byte NO_BASE = 0;

    private void buildSequence()
    {
        int currentIndex = mIsForward ? mBases.length - 1 : 0;
        int refPosition = mJunctionPosition;

        // set the junction read base as established from the initial split reads
        mBases[currentIndex] = mAssembly.bases()[mAssembly.junctionIndex()];
        mBaseQuals[currentIndex] = mAssembly.baseQuals()[mAssembly.junctionIndex()];

        int baseCount = Nucleotides.DNA_BASES.length;

        CigarOperator currentElementType = M;
        int currentElementLength = 0;

        List<RefReadParseState> activeReads = mReads.stream().filter(x -> x.isValid()).collect(Collectors.toList());

        // boolean isIndelJunction = mAssembly.junction().indelBased();

        while(!activeReads.isEmpty())
        {
            if(currentElementType == M || currentElementType == I)
                currentIndex += mIsForward ? -1 : 1;

            ++currentElementLength;

            if(currentIndex < 0 || currentIndex >= mBases.length)
                break;

            byte consensusBase = 0;
            int consensusMaxQual = 0;
            int consensusQualTotal = 0;
            int consensusReadCount = 0;

            // per-base arrays are only used for high-qual mismatches
            int[] readCounts = null;
            int[] totalQuals = null;
            int[] maxQuals = null;

            // move to the next position or index if during an insert
            progressReadState(activeReads, currentElementType);

            if(activeReads.isEmpty())
                break;

            // establish new most common operator
            CigarOperator nextElementType = findNextOperator(activeReads);

            if(nextElementType != currentElementType)
            {
                mCigarElements.add(new CigarElement(currentElementLength, currentElementType));
                currentElementLength = 0;
                currentElementType = nextElementType;
            }

            if(currentElementType == D)
            {
                refPosition += mIsForward ? -1 : 1;
                markReadBaseMatches(activeReads, currentElementType, NO_BASE, NO_BASE);
                continue;
            }

            // now establish the consensus base
            for(RefReadParseState read : activeReads)
            {
                if(read.operator() != currentElementType)
                    continue;

                byte base = read.currentBase();
                int qual = read.currentQual();

                if(base == DNA_N_BYTE)
                {
                    base = DNA_BASE_BYTES[0];
                    qual = 0;
                }

                if(readCounts == null)
                {
                    if(consensusBase == NO_BASE || (base != consensusBase && belowMinQual(consensusMaxQual) && aboveMinQual(qual)))
                    {
                        // set first or replace with first high qual
                        consensusBase = base;
                        consensusMaxQual = qual;
                        consensusQualTotal = qual;
                        consensusReadCount = 1;
                        continue;
                    }
                    else if(base == consensusBase)
                    {
                        consensusMaxQual = max(qual, consensusMaxQual);
                        consensusQualTotal += qual;
                        ++consensusReadCount;
                        continue;
                    }
                    else if(base != consensusBase && belowMinQual(qual))
                    {
                        // low-qual disagreement - ignore regardless of consensus qual
                        continue;
                    }

                    // high-qual mismatch so start tracking frequencies for each base
                    readCounts = new int[baseCount];
                    totalQuals = new int[baseCount];
                    maxQuals = new int[baseCount];

                    // back port existing counts to the per-base arrays
                    int baseIndex = Nucleotides.baseIndex(consensusBase);
                    totalQuals[baseIndex] = consensusQualTotal;
                    maxQuals[baseIndex] = consensusMaxQual;
                    readCounts[baseIndex] = consensusReadCount;
                }

                int baseIndex = Nucleotides.baseIndex(base);

                totalQuals[baseIndex] += qual;
                maxQuals[baseIndex] = max(maxQuals[baseIndex], qual);
                ++readCounts[baseIndex];
            }

            if(readCounts != null)
            {
                // take the bases with the highest qual totals
                int maxQual = 0;
                int maxBaseIndex = 0;
                for(int b = 0; b < baseCount; ++b)
                {
                    if(totalQuals[b] > maxQual)
                    {
                        maxQual = totalQuals[b];
                        maxBaseIndex = b;
                    }
                }

                consensusBase = DNA_BASE_BYTES[maxBaseIndex];
                consensusMaxQual = (byte)maxQuals[maxBaseIndex];
            }

            mBases[currentIndex] = consensusBase;
            mBaseQuals[currentIndex] = (byte)consensusMaxQual;

            // mark active reads as matching or not
            markReadBaseMatches(activeReads, currentElementType, consensusBase, consensusMaxQual);

            if(currentElementType == M || currentElementType == D)
                refPosition += mIsForward ? -1 : 1;
        }

        // add the last, current element
        mCigarElements.add(new CigarElement(currentElementLength, currentElementType));

        mRefBasePosition = refPosition;

        // correct ref position where M bases have been skipped over from an insert
        if(mCigarElements.stream().anyMatch(x -> x.getOperator() == I && x.getLength() >= MIN_INDEL_SUPPORT_LENGTH))
        {
            int outerRefPosition = -1;

            for(RefReadParseState read : mReads)
            {
                if(!read.isValid())
                    continue;

                if(mIsForward)
                {
                    if(outerRefPosition < 0 || read.refPosition() < outerRefPosition)
                        outerRefPosition = read.refPosition();
                }
                else
                {
                    outerRefPosition = max(outerRefPosition, read.refPosition());
                }
            }

            mRefBasePosition = outerRefPosition;
        }

        if(mIsForward)
            Collections.reverse(mCigarElements);
    }

    private static CigarOperator findNextOperator(final List<RefReadParseState> reads)
    {
        int inserts = 0;
        int deletes = 0;
        int aligned = 0;

        for(RefReadParseState read : reads)
        {
            if(read.operator() == M)
                ++aligned;
            else if(read.operator() == D)
                ++deletes;
            else if(read.operator() == I)
                ++inserts;
        }

        if(aligned >= inserts && aligned >= deletes)
            return M;

        return inserts >= deletes ? I : D;
    }

    private static void progressReadState(final List<RefReadParseState> reads, final CigarOperator currentElementType)
    {
        // move to the next position or index if during an insert
        boolean allowNonIndelMatch = false;

        if(currentElementType == I)
        {
            int indelCount = 0;

            for(RefReadParseState read : reads)
            {
                if(read.operator() == I && read.elementLength() >= MIN_INDEL_SUPPORT_LENGTH)
                {
                    ++indelCount;

                    if(indelCount > reads.size() / 2)
                    {
                        allowNonIndelMatch = true;
                        break;
                    }
                }
            }
        }

        int index = 0;
        while(index < reads.size())
        {
            RefReadParseState read = reads.get(index);

            if(currentElementType == M || currentElementType == D)
            {
                if(read.operator() == I && read.elementLength() < MIN_INDEL_SUPPORT_LENGTH)
                {
                    read.skipInsert();
                }
                else
                {
                    read.moveNext();
                }
            }
            else // insert
            {
                if(read.operator() == I || (read.operator() == M && allowNonIndelMatch))
                {
                    read.moveNext();
                }
            }

            if(read.exhausted() || read.exceedsMaxMismatches())
                reads.remove(index);
            else
                ++index;
        }
    }

    private static void markReadBaseMatches(
            final List<RefReadParseState> reads, final CigarOperator currentElementType, final byte consensusBase, final int consensusQual)
    {
        for(RefReadParseState read : reads)
        {
            if(currentElementType == D)
            {
                if(read.operator() != currentElementType)
                    read.addIndelMismatch();

                continue;
            }

            boolean aboveMinQual = aboveMinQual(read.currentQual());
            boolean consensusAboveMinQual = aboveMinQual(consensusQual);
            boolean operatorMismatch = read.operator() != currentElementType;

            if(operatorMismatch && read.operator() == D)
            {
                read.addIndelMismatch();
                continue;
            }

            if(read.currentBase() != consensusBase)
            {
                if(aboveMinQual && consensusAboveMinQual)
                {
                    if(operatorMismatch)
                        read.addIndelMismatch();
                    else
                        read.addMismatch();
                }
            }
            else
            {
                if(aboveMinQual)
                    read.addHighQualMatch();
            }
        }
    }

    private void trimFinalSequence()
    {
        int currentIndex = mIsForward ? mBases.length - 1 : 0;
        int validLength = 0;

        while(currentIndex >= 0 && currentIndex < mBases.length)
        {
            if(mBases[currentIndex] == 0)
                break;

            ++validLength;
            currentIndex += mIsForward ? -1 : 1;
        }

        if(validLength == mBases.length)
            return;

        int reduction = mBases.length - validLength;
        int startIndex = mIsForward ? reduction : 0;
        int endIndex = startIndex + validLength - 1;
        mBases = subsetArray(mBases, startIndex, endIndex);
        mBaseQuals = subsetArray(mBaseQuals, startIndex, endIndex);
    }

    public String toString()
    {
        return format("refPosition(%d len=%d) cigar(%s) reads(%s)",
                mRefBasePosition, mBases.length, cigarStr(), mReads.size());
    }

    @VisibleForTesting
    public String refBaseSequence() { return new String(mBases); }
}