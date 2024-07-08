package com.hartwig.hmftools.sage.testutil;

import static java.lang.Math.floor;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.google.common.primitives.Ints.min;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.X;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class MutatedBases
{
    public static class MutatedBase
    {
        public final int RefPos;
        public final char Base;
        public final CigarOperator CigarOp;

        public MutatedBase(int refPos, char base, final CigarOperator cigarOp)
        {
            RefPos = refPos;
            Base = base;
            CigarOp = cigarOp;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof MutatedBase))
            {
                return false;
            }
            final MutatedBase that = (MutatedBase) o;
            return RefPos == that.RefPos && Base == that.Base && CigarOp == that.CigarOp;
        }
    }

    @VisibleForTesting
    public static class AlignedBases
    {
        public final int LeftSoftclipLength;
        public final int RightSoftclipLength;
        public final List<MutatedBase> Aligned;

        public AlignedBases(final int leftSoftclipLength, final int rightSoftclipLength, final List<MutatedBase> aligned)
        {
            LeftSoftclipLength = leftSoftclipLength;
            RightSoftclipLength = rightSoftclipLength;
            Aligned = aligned;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof AlignedBases))
            {
                return false;
            }
            final AlignedBases that = (AlignedBases) o;
            return LeftSoftclipLength == that.LeftSoftclipLength && RightSoftclipLength == that.RightSoftclipLength
                    && Objects.equals(Aligned, that.Aligned);
        }
    }

    private final String mRefBases;
    private final List<MutatedBase> mMutatedBases;
    private final BaseRegion mMutationBounds;

    public MutatedBases(final String refBases, final List<MutatedBase> mutatedBases, @Nullable final BaseRegion mutationBounds)
    {
        mRefBases = refBases;
        mMutatedBases = mutatedBases;
        mMutationBounds = mutationBounds;
    }

    @VisibleForTesting
    public MutatedBases(final String refBases, final List<MutatedBase> mutatedBases)
    {
        this(refBases, mutatedBases, null);
    }

    @Nullable
    public SAMRecord getPairedRead(final String readName, final String chromosome, int mutStartIndex, int readLength, int fragmentLength,
            boolean firstInPair, boolean onForwardStrand)
    {
        List<MutatedBase> readBases = mMutatedBases.subList(mutStartIndex, mutStartIndex + readLength);
        AlignedBases aligned = getAlignedBases(mRefBases, readBases);
        if(aligned == null)
        {
            return null;
        }

        String readBasesStr = readBases.stream().map(x -> String.valueOf(x.Base)).collect(Collectors.joining());
        int alignmentStart = aligned.Aligned.get(0).RefPos;
        String cigarStr = getCigarStr(aligned);
        String mateCigarStr = String.valueOf(readLength) + "M";
        int mateAlignmentStart;
        if(onForwardStrand)
        {
            mateAlignmentStart = mMutatedBases.get(mutStartIndex + fragmentLength - readLength).RefPos;
        }
        else
        {
            mateAlignmentStart = mMutatedBases.get(mutStartIndex + readLength - fragmentLength).RefPos;
        }

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readName, chromosome, alignmentStart, readBasesStr, cigarStr, chromosome, mateAlignmentStart, !onForwardStrand,
                false, null, onForwardStrand, mateCigarStr);

        read.setFirstOfPairFlag(firstInPair);
        read.setSecondOfPairFlag(!firstInPair);

        return read;
    }

    @Nullable
    public SAMRecord getUnpairedRead(final String readName, final String chromosome, int mutStartIndex, int readLength, boolean onForwardStrand)
    {
        List<MutatedBase> readBases = mMutatedBases.subList(mutStartIndex, mutStartIndex + readLength);
        AlignedBases aligned = getAlignedBases(mRefBases, readBases);
        if(aligned == null)
        {
            return null;
        }

        String readBasesStr = readBases.stream().map(x -> String.valueOf(x.Base)).collect(Collectors.joining());
        int alignmentStart = aligned.Aligned.get(0).RefPos;
        String cigarStr = getCigarStr(aligned);

        SAMRecord read = createSamRecordUnpaired(readName, chromosome, alignmentStart, readBasesStr, cigarStr, !onForwardStrand, false, null);

        return read;
    }

    private static final int DEFAULT_MAX_FAILURES = 1000;

    public List<SAMRecord> generateRandomPairedReads(final Random random, String chromosome, int depth, int readLength, int fragmentLength,
            final BaseRegion overlappingRegion)
    {
        return generateRandomPairedReads(random, chromosome, depth, readLength, fragmentLength, overlappingRegion, DEFAULT_MAX_FAILURES);
    }

    private static double averageOverlappingBases(int readLength, int regionLength)
    {
        int overlappingBasesCount = 0;
        for(int i = 1; i <= regionLength; ++i)
        {
            overlappingBasesCount += min(i, regionLength - i + 1, readLength);
        }

        return 1.0 * overlappingBasesCount / regionLength;
    }

    public List<SAMRecord> generateRandomPairedReads(final Random random, String chromosome, int depth, int readLength, int fragmentLength,
            final BaseRegion region, int maxFailures)
    {
        List<SAMRecord> reads = Lists.newArrayList();
        int regionStartMutIndex = region.start() - readLength + 1;
        int regionEndMutIndex = region.end() + readLength - 1;
        int totalBases = depth * (regionEndMutIndex - regionStartMutIndex + 1);
        double avgOverlappingBases = averageOverlappingBases(readLength, region.baseLength());
        int totalReads = (int) floor(round(totalBases / avgOverlappingBases));

        int failureCount = 0;
        int totalReadsDigitCount = String.valueOf(totalReads).length();
        String readNameFormat = format("READ_%%0%dd", totalReadsDigitCount);
        for(int i = 0; i < totalReads && failureCount <= maxFailures; )
        {
            String readName = format(readNameFormat, i);
            int mutPosStart = regionStartMutIndex + random.nextInt(regionEndMutIndex - regionStartMutIndex - readLength + 1);
            boolean firstInPair = random.nextBoolean();
            boolean onForwardStrand = random.nextBoolean();
            SAMRecord read = getPairedRead(readName, chromosome, mutPosStart, readLength, fragmentLength, firstInPair, onForwardStrand);
            if(read == null)
            {
                ++failureCount;
                continue;
            }

            reads.add(read);
            ++i;
        }

        if(failureCount > maxFailures)
        {
            throw new RuntimeException("Max failures exceeds when generating random reads");
        }

        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        return reads;
    }

    public int leftMutIndexFromRefPos(int refPos)
    {
        if(mMutatedBases.isEmpty())
        {
            return -1;
        }

        MutatedBase key = new MutatedBase(refPos, '.', M);
        int bisectIndex = Collections.binarySearch(mMutatedBases, key, Comparator.comparingInt(a -> a.RefPos));

        if(bisectIndex >= 0)
        {
            // exact match.
            return bisectIndex;
        }

        int insertionPoint = -bisectIndex - 1;
        return insertionPoint - 1;
    }

    public int rightMutIndexFromRefPos(int refPos)
    {
        if(mMutatedBases.isEmpty())
        {
            return -1;
        }

        MutatedBase key = new MutatedBase(refPos, '.', M);
        int bisectIndex = Collections.binarySearch(mMutatedBases, key, Comparator.comparingInt(a -> a.RefPos));

        if(bisectIndex >= 0)
        {
            // exact match.
            return bisectIndex;
        }

        int insertionPoint = -bisectIndex - 1;
        if(insertionPoint >= mMutatedBases.size())
        {
            return -1;
        }

        return insertionPoint;
    }

    @Nullable
    public BaseRegion mutIndexMutationBounds()
    {
        if(mMutationBounds == null)
            return null;

        int startMutIndex = leftMutIndexFromRefPos(mMutationBounds.start());
        if(startMutIndex < 0)
            return null;

        int endMutIndex = rightMutIndexFromRefPos(mMutationBounds.end());
        if(endMutIndex < 0)
            return null;

        return new BaseRegion(startMutIndex, endMutIndex);
    }

    @VisibleForTesting
    @Nullable
    public BaseRegion refPosMutationBounds()
    {
        return mMutationBounds;
    }

    @VisibleForTesting
    @Nullable
    public static AlignedBases getAlignedBases(final String refBases, final List<MutatedBase> readBases)
    {
        if(readBases.isEmpty())
        {
            return null;
        }

        Deque<MutatedBase> leftScBases = new ArrayDeque<>();
        Deque<MutatedBase> alignmentBases = new ArrayDeque<>(readBases);
        Deque<MutatedBase> rightScBases = new ArrayDeque<>();

        // initialise splits by first and last matches
        while(!alignmentBases.isEmpty() && alignmentBases.peekFirst().CigarOp != EQ)
        {
            leftScBases.addLast(alignmentBases.pollFirst());
        }

        while(!alignmentBases.isEmpty() && alignmentBases.peekLast().CigarOp != EQ)
        {
            rightScBases.addFirst(alignmentBases.pollLast());
        }

        if(alignmentBases.isEmpty())
        {
            return null;
        }

        // absorb extra matches into the alignment
        while(!leftScBases.isEmpty())
        {
            MutatedBase base = leftScBases.peekLast();
            MutatedBase nextBase = alignmentBases.peekFirst();
            if(base.Base != refBases.charAt(nextBase.RefPos - 2))
            {
                break;
            }

            leftScBases.removeLast();
            alignmentBases.addFirst(new MutatedBase(nextBase.RefPos - 1, base.Base, EQ));
        }

        while(!rightScBases.isEmpty())
        {
            MutatedBase base = rightScBases.peekFirst();
            MutatedBase prevBase = alignmentBases.peekLast();
            if(base.Base != refBases.charAt(prevBase.RefPos))
            {
                break;
            }

            rightScBases.removeFirst();
            alignmentBases.addLast(new MutatedBase(prevBase.RefPos + 1, base.Base, EQ));
        }

        return new AlignedBases(leftScBases.size(), rightScBases.size(), Collections.unmodifiableList(Lists.newArrayList(alignmentBases)));
    }

    @VisibleForTesting
    public static String getCigarStr(final List<MutatedBase> bases)
    {
        if(bases.isEmpty())
        {
            return "";
        }

        UnaryOperator<CigarOperator> cigarOpMapper = (final CigarOperator op) -> op == X || op == EQ ? M : op;

        StringBuilder cigarBuilder = new StringBuilder();
        MutatedBase prevBase = bases.get(0);
        CigarOperator currentOp = cigarOpMapper.apply(prevBase.CigarOp);
        int currentLength = 1;
        for(int i = 1; i < bases.size(); ++i)
        {
            MutatedBase base = bases.get(i);
            CigarOperator op = cigarOpMapper.apply(base.CigarOp);
            int delLength = base.RefPos - prevBase.RefPos - 1;
            if(delLength > 0)
            {
                if(currentLength > 0)
                {
                    cigarBuilder.append(currentLength);
                    cigarBuilder.append(currentOp.name());
                }

                cigarBuilder.append(delLength);
                cigarBuilder.append(D.name());

                currentOp = null;
                currentLength = 0;
            }

            if(currentLength == 0)
            {
                currentOp = op;
                currentLength = 1;
            }
            else if(op == currentOp)
            {
                ++currentLength;
            }
            else
            {
                cigarBuilder.append(currentLength);
                cigarBuilder.append(currentOp.name());

                currentOp = op;
                currentLength = 1;
            }

            prevBase = base;
        }

        if(currentLength > 0)
        {
            cigarBuilder.append(currentLength);
            cigarBuilder.append(currentOp.name());
        }

        return cigarBuilder.toString();
    }

    @VisibleForTesting
    public static String getCigarStr(final AlignedBases alignedBases)
    {
        StringBuilder cigarBuilder = new StringBuilder();
        if(alignedBases.LeftSoftclipLength > 0)
        {
            cigarBuilder.append(alignedBases.LeftSoftclipLength);
            cigarBuilder.append('S');
        }

        cigarBuilder.append(getCigarStr(alignedBases.Aligned));

        if(alignedBases.RightSoftclipLength > 0)
        {
            cigarBuilder.append(alignedBases.RightSoftclipLength);
            cigarBuilder.append('S');
        }

        return cigarBuilder.toString();
    }

    @VisibleForTesting
    public List<MutatedBase> mutatedBases()
    {
        return mMutatedBases;
    }
}
