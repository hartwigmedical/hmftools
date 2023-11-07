package com.hartwig.hmftools.markdups.consensus;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.SAMRecord;

public class BaseBuilder
{
    private final RefGenomeInterface mRefGenome;

    // cached for the majority of successive reads being on the same chromosome, to protect against ref genome base requests beyond limits
    private int mChromosomeLength;

    public BaseBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLength = 0;
    }

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

        for(int i = 0; i < readCount; ++i)
        {
            readOffsets[i] = reads.get(i).getReadBases().length - baseLength;
        }


        int firstInPairCount = (int) reads.stream().filter(q -> q.getFirstOfPairFlag()).count();
        int secondInPairCount = readCount - firstInPairCount;
        boolean isDual = firstInPairCount > 0 && secondInPairCount > 0;
        boolean[] isFirstInPair = new boolean[readCount];;
        for(int i = 0; i < readCount; ++i)
            isFirstInPair[i] = reads.get(i).getFirstOfPairFlag();

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];
        // we also track bases and quals by first/second in pair status
        byte[][] locationBasesPair = new byte[][] { new byte[readCount], new byte[readCount] };
        byte[][] locationQualsPair = new byte[][] { new byte[readCount], new byte[readCount] };

        boolean dualMismatchedLogged = false;

        for(int baseIndex = 0; baseIndex < baseLength; ++baseIndex)
        {
            // check bases at this index
            // work on the premise that most bases will agree
            boolean hasMismatch = false;
            int maxQual = 0;
            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                SAMRecord read = reads.get(r);

                locationBases[r] = NO_BASE;
                locationQuals[r] = 0;
                locationBasesPair[0][r] = NO_BASE;
                locationBasesPair[1][r] = NO_BASE;
                locationQualsPair[0][r] = 0;
                locationQualsPair[1][r] = 0;


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

                int pairIdx = isFirstInPair[r] ? 0 : 1;
                locationBasesPair[pairIdx][r] = locationBases[r];
                locationQualsPair[pairIdx][r] = locationQuals[r];

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

                if (!dualMismatchedLogged && isDual)
                {
                    MD_LOGGER.trace("Mismatched basis found in dual stranded read group readCount({})", reads.size());
                    for (SAMRecord read : reads)
                        MD_LOGGER.trace("Read in dual stranded read group: {}", read);

                    dualMismatchedLogged = true;
                }

                byte[] consensusBaseAndQual;
                if(isDual && basePosition != -1)
                {
                    consensusBaseAndQual =
                            determineBaseAndQualDual(locationBasesPair, locationQualsPair, chromosome, basePosition);
                }
                else
                {
                    consensusBaseAndQual = determineBaseAndQual(locationBases, locationQuals, chromosome, basePosition);
                }

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }
        }
    }

    public byte[] determineBaseAndQualDual(
            final byte[][] locationBasesPair, final byte[][] locationQualsPair, final String chromosome, int position)
    {
        byte[] firstConsensusBaseAndQual = determineBaseAndQual(locationBasesPair[0], locationQualsPair[0], chromosome, position);
        byte[] secondConsensusBaseAndQual = determineBaseAndQual(locationBasesPair[1], locationQualsPair[1], chromosome, position);

        if(firstConsensusBaseAndQual[0] == NO_BASE)
        {
            return secondConsensusBaseAndQual;
        }

        if(secondConsensusBaseAndQual[0] == NO_BASE)
        {
            return firstConsensusBaseAndQual;
        }

        if(firstConsensusBaseAndQual[0] == secondConsensusBaseAndQual[0])
        {
            byte qual = (byte) max(firstConsensusBaseAndQual[1], secondConsensusBaseAndQual[1]);
            return new byte[] { firstConsensusBaseAndQual[0], qual };
        }

        byte refBase = mRefGenome.getBaseString(chromosome, position, position).getBytes()[0];
        boolean firstIsRef = firstConsensusBaseAndQual[0] == refBase;
        boolean secondIsRef = secondConsensusBaseAndQual[0] == refBase;
        if(!firstIsRef && !secondIsRef)
        {
            byte[] locationBasePair = new byte[] { firstConsensusBaseAndQual[0], secondConsensusBaseAndQual[0] };
            byte[] locationQualPair = new byte[] { firstConsensusBaseAndQual[1], secondConsensusBaseAndQual[1] };
            return determineBaseAndQual(locationBasePair, locationQualPair, chromosome, position);
        }

        int refQual;
        int differingQual;
        if(firstIsRef)
        {
            refQual = firstConsensusBaseAndQual[1];
            differingQual = secondConsensusBaseAndQual[1];
        }
        else
        {
            refQual = secondConsensusBaseAndQual[1];
            differingQual = firstConsensusBaseAndQual[1];
        }

        byte qual = (byte) max(0.0, refQual - differingQual);
        return new byte[] { refBase, qual };
    }

    public byte[] determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position)
    {
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

        int differingQual = 0;

        for(int i = 0; i < distinctBases.size(); ++i)
        {
            if(distinctBases.get(i) != maxBase)
                differingQual += qualTotals.get(i);
        }

        double calcQual = (double)maxQual * max(0.0, maxQualTotal - differingQual) / maxQualTotal;

        return new byte[] { maxBase, (byte)round(calcQual) };
    }

    public void setChromosomLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }
    public int chromosomeLength() { return mChromosomeLength; }

    public RefGenomeInterface refGenome() { return mRefGenome; }
}
