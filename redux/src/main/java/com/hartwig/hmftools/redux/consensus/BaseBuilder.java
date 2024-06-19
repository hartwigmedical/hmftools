package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.BASE_QUAL_MINIMUM;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.qual.BaseQualAdjustment;

import htsjdk.samtools.SAMRecord;

public class BaseBuilder
{
    private final RefGenomeInterface mRefGenome;
    private final ConsensusStatistics mConsensusStats;

    // cached for the majority of successive reads being on the same chromosome, to protect against ref genome base requests beyond limits
    private int mChromosomeLength;

    public BaseBuilder(final RefGenomeInterface refGenome, final ConsensusStatistics consensusStats)
    {
        mRefGenome = refGenome;
        mConsensusStats = consensusStats;
        mChromosomeLength = 0;
    }

    public void setChromosomLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }
    public int chromosomeLength() { return mChromosomeLength; }
    public RefGenomeInterface refGenome() { return mRefGenome; }

    public static final byte NO_BASE = 0;
    public static final int INVALID_POSITION = -1;

    public void buildReadBases(final List<SAMRecord> reads, final ConsensusState consensusState)
    {
        int chromosomeLength = mChromosomeLength;
        if(chromosomeLength == 0)
            chromosomeLength = mRefGenome.getChromosomeLength(reads.get(0).getReferenceName());

        int baseLength = consensusState.Bases.length;
        int readCount = reads.size();
        String chromosome = reads.get(0).getContig();

        int[] readOffsets = new int[readCount];
        boolean[] isFirstInPair = new boolean[readCount];
        boolean isDualStrand = isDualStrandAndIsFirstInPair(reads, isFirstInPair);

        for(int i = 0; i < readCount; ++i)
        {
            readOffsets[i] = reads.get(i).getReadBases().length - baseLength;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        for(int baseIndex = 0; baseIndex < baseLength; ++baseIndex)
        {
            // check bases at this index - work on the premise that most bases will agree
            boolean hasMismatch = false;
            int maxQual = 0;
            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                SAMRecord read = reads.get(r);

                locationBases[r] = NO_BASE;

                int readIndex;
                if(consensusState.IsForward)
                {
                    readIndex = baseIndex;

                    if(readOffsets[r] != 0 && baseIndex >= read.getReadBases().length)
                        continue;
                }
                else
                {
                    readIndex = baseIndex + readOffsets[r];

                    if(readIndex < 0)
                        continue;
                }

                locationBases[r] = reads.get(r).getReadBases()[readIndex];
                locationQuals[r] = reads.get(r).getBaseQualities()[readIndex];

                if(firstBase == NO_BASE)
                    firstBase = locationBases[r];
                else
                    hasMismatch |= locationBases[r] != firstBase;

                maxQual = max(locationQuals[r], maxQual);
            }

            if(!hasMismatch)
            {
                consensusState.Bases[baseIndex] = firstBase;
                consensusState.BaseQualities[baseIndex] = (byte)maxQual;
            }
            else
            {
                int basePosition = consensusState.MinUnclippedPosStart + baseIndex;

                if(basePosition < 1 || basePosition > chromosomeLength)
                    basePosition = INVALID_POSITION; // protect against over-runs from soft-clips - rare but possible

                byte[] consensusBaseAndQual;

                if(isDualStrand && basePosition != INVALID_POSITION)
                {
                    // split the reads into 2 consensus reads and then compare
                    consensusBaseAndQual = determineDualStrandBaseAndQual(isFirstInPair, locationBases, locationQuals, chromosome, basePosition);
                }
                else
                {
                    consensusBaseAndQual = determineBaseAndQual(locationBases, locationQuals, chromosome, basePosition);
                }

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = BaseQualAdjustment.adjustBaseQual(consensusBaseAndQual[1]);
            }
        }
    }

    public byte[] determineDualStrandBaseAndQual(
            final boolean[] isFirstInPair, final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position)
    {
        int readCount = isFirstInPair.length;

        int firstInPairCount = 0;
        int secondInPairCount = 0;

        for(int i = 0; i < isFirstInPair.length; ++i)
        {
            if(isFirstInPair[i])
                ++firstInPairCount;
            else
                ++secondInPairCount;
        }

        byte[] locationBasesFirst = new byte[firstInPairCount];
        byte[] locationQualsFirst = new byte[firstInPairCount];

        byte[] locationBasesSecond = new byte[secondInPairCount];
        byte[] locationQualsSecond = new byte[secondInPairCount];

        int firstIndex = 0;
        int secondIndex = 0;

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(isFirstInPair[i])
            {
                locationBasesFirst[firstIndex] = locationBases[i];
                locationQualsFirst[firstIndex] = locationQuals[i];
                ++firstIndex;
            }
            else
            {
                locationBasesSecond[secondIndex] = locationBases[i];
                locationQualsSecond[secondIndex] = locationQuals[i];
                ++secondIndex;
            }
        }

        byte[] firstBaseAndQual = determineBaseAndQual(locationBasesFirst, locationQualsFirst, chromosome, position);
        byte[] secondBaseAndQual = determineBaseAndQual(locationBasesSecond, locationQualsSecond, chromosome, position);

        if(firstBaseAndQual[0] == NO_BASE)
            return secondBaseAndQual;

        if(secondBaseAndQual[0] == NO_BASE)
            return firstBaseAndQual;

        if(firstBaseAndQual[0] == secondBaseAndQual[0])
        {
            byte qual = (byte)max(firstBaseAndQual[1], secondBaseAndQual[1]);
            return new byte[] { firstBaseAndQual[0], qual };
        }

        mConsensusStats.registerDualStrandMismatchReadGroup(readCount);

        byte refBase = mRefGenome.getBaseString(chromosome, position, position).getBytes()[0];
        boolean firstIsRef = firstBaseAndQual[0] == refBase;
        boolean secondIsRef = secondBaseAndQual[0] == refBase;

        if(!firstIsRef && !secondIsRef)
        {
            byte maxBase;
            int maxQual;
            int differingQual;
            if(firstBaseAndQual[1] >= secondBaseAndQual[1])
            {
                maxBase = firstBaseAndQual[0];
                maxQual = firstBaseAndQual[1];
                differingQual = secondBaseAndQual[1];
            }
            else
            {
                maxBase = secondBaseAndQual[0];
                maxQual = secondBaseAndQual[1];
                differingQual = firstBaseAndQual[1];
            }

            byte qual = (byte) max(0, maxQual - differingQual);
            return new byte[] { maxBase, qual };
        }

        int refQual;
        int differingQual;
        if(firstIsRef)
        {
            refQual = firstBaseAndQual[1];
            differingQual = secondBaseAndQual[1];
        }
        else
        {
            refQual = secondBaseAndQual[1];
            differingQual = firstBaseAndQual[1];
        }

        byte qual = (byte) max(BASE_QUAL_MINIMUM, refQual - differingQual);
        return new byte[] { refBase, qual };
    }

    public byte[] determineBaseAndQual(final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position)
    {
        if(locationBases.length == 1)
        {
            // early exit for dual strand with a single read on one side - a very common scenario
            return new byte[] { locationBases[0], locationQuals[0] };
        }

        List<Byte> distinctBases = Lists.newArrayListWithCapacity(4);
        List<Integer> qualTotals = Lists.newArrayListWithCapacity(4);
        List<Integer> maxQuals = Lists.newArrayListWithCapacity(4);

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == NO_BASE)
                continue;

            boolean found = false;

            for(int j = 0; j < distinctBases.size(); ++j)
            {
                if(distinctBases.get(j) == locationBases[i])
                {
                    int qualTotal = qualTotals.get(j) + locationQuals[i];
                    qualTotals.set(j, qualTotal);
                    maxQuals.set(j, max(maxQuals.get(j), locationQuals[i]));
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                distinctBases.add(locationBases[i]);
                qualTotals.add((int)locationQuals[i]);
                maxQuals.add((int)locationQuals[i]);
            }
        }

        if(distinctBases.isEmpty())
        {
            return new byte[] { NO_BASE, 0 };
        }

        byte maxBase = distinctBases.get(0);
        boolean maxIsRef = false;
        int maxQual = maxQuals.get(0);
        int maxQualTotal = qualTotals.get(0);

        for(int i = 1; i < distinctBases.size(); ++i)
        {
            if(qualTotals.get(i) > maxQualTotal)
            {
                maxQualTotal = qualTotals.get(i);
                maxQual = maxQuals.get(i);
                maxBase = distinctBases.get(i);
            }
            else if(chromosome != null && qualTotals.get(i) >= maxQualTotal && !maxIsRef && position != INVALID_POSITION)
            {
                // chromosome will be null for unmapped reads
                String refBase = mRefGenome.getBaseString(chromosome, position, position);

                if(maxBase == refBase.getBytes()[0])
                {
                    maxIsRef = true;
                }
                else if(distinctBases.get(i) == refBase.getBytes()[0])
                {
                    maxQualTotal = qualTotals.get(i);
                    maxQual = maxQuals.get(i);
                    maxBase = distinctBases.get(i);
                    maxIsRef = true;
                }
            }
        }

        // collect base quals matching the selected base to find the median
        List<Integer> selectBaseQuals = Lists.newArrayList();

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == maxBase)
                selectBaseQuals.add((int)locationQuals[i]);
        }

        Collections.sort(selectBaseQuals);

        int medianBaseQualIndex = selectBaseQuals.size() / 2;
        int medianBaseQual = selectBaseQuals.get(medianBaseQualIndex);

        int differingQual = 0;

        for(int i = 0; i < distinctBases.size(); ++i)
        {
            if(distinctBases.get(i) != maxBase)
                differingQual += qualTotals.get(i);
        }

        double calcQual = (double)medianBaseQual * max(BASE_QUAL_MINIMUM, maxQualTotal - differingQual) / maxQualTotal;

        return new byte[] { maxBase, (byte)round(calcQual) };
    }

    public static boolean isDualStrandAndIsFirstInPair(final List<SAMRecord> reads, final boolean[] isFirstInPairOut)
    {
        boolean isDualStrand = false;

        for(int i = 0; i < reads.size(); ++i)
        {
            isFirstInPairOut[i] = reads.get(i).getFirstOfPairFlag();

            if(i > 0)
                isDualStrand |= isFirstInPairOut[0] != isFirstInPairOut[i];
        }

        return isDualStrand;
    }
}
