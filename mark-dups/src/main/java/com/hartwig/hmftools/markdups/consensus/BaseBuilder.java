package com.hartwig.hmftools.markdups.consensus;

import static java.lang.Math.max;
import static java.lang.Math.round;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.SAMRecord;

public class BaseBuilder
{
    private final RefGenomeInterface mRefGenome;

    private int mChromosomeLength;

    public BaseBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLength = 0;
    }

    public static final byte NO_BASE = 0;

    public void buildReadBases(final List<SAMRecord> reads, final ConsensusState consensusState)
    {
        int chromosomeLength = mChromosomeLength;
        if(chromosomeLength == 0)
            chromosomeLength = mRefGenome.getChromosomeLength(reads.get(0).getReferenceName());

        int baseLength = consensusState.Bases.length;

        int readCount = reads.size();

        int[] readOffsets = new int[readCount];

        for(int i = 0; i < readCount; ++i)
        {
            readOffsets[i] = reads.get(i).getReadBases().length - baseLength;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

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

                locationBases[r] = 0;

                int readIndex;
                if(consensusState.IsForward)
                {
                    readIndex = baseIndex;

                    if(readOffsets[r] != 0 && baseIndex >= read.getReadBases().length)
                    {
                        locationBases[r] = NO_BASE;
                        locationQuals[r] = NO_BASE;
                        continue;
                    }
                }
                else
                {
                    readIndex = baseIndex + readOffsets[r];

                    if(readIndex < 0)
                    {
                        locationBases[r] = NO_BASE;
                        locationQuals[r] = NO_BASE;
                        continue;
                    }
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
                if(basePosition > chromosomeLength)
                    basePosition = -1;

                byte[] consensusBaseAndQual = determineBaseAndQual(
                        locationBases, locationQuals, reads.get(0).getContig(), basePosition);

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }
        }
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
            else if(chromosome != null && qualTotals.get(i) >= maxQualTotal && !maxIsRef
                    && position != -1) // chromosome will be null for unmapped reads
            {
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

    public void setChromosomLength(int chromosomeLength)
    {
        mChromosomeLength = chromosomeLength;
    }

    public int chromosomeLength()
    {
        return mChromosomeLength;
    }

    public RefGenomeInterface refGenome()
    {
        return mRefGenome;
    }
}
