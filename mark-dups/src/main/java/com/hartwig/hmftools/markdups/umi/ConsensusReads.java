package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.UNSET;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.SAMRecord;

public class ConsensusReads
{
    private final UmiConfig mConfig;
    private final RefGenomeInterface mRefGenome;

    public ConsensusReads(final UmiConfig config, final RefGenomeInterface refGenome)
    {
        mConfig = config;
        mRefGenome = refGenome;
    }

    public ConsensusReadInfo createConsensusRead(final List<SAMRecord> reads, final String groupIdentifier)
    {
        if(reads.size() <= 1)
            return null;

        boolean isForward = !reads.get(0).getReadNegativeStrandFlag();
        int maxBaseLength = 0;
        boolean hasIndels = false;

        // work out the outermost boundaries - soft-clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState();

        for(SAMRecord read : reads)
        {
            maxBaseLength = max(maxBaseLength, read.getReadBases().length);

            hasIndels |= read.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == I || x.getOperator() == D);

            consensusState.setBoundaries(read);
        }

        consensusState.setBaseLength(maxBaseLength);

        if(hasIndels)
        {

        }
        else
        {
            buildReadComponents(reads, isForward, consensusState);

        }

        SAMRecord consensusRead = consensusState.createConsensusRead(reads.get(0), groupIdentifier);

        return new ConsensusReadInfo(consensusRead, consensusState.outcome());
    }

    private void buildReadComponents(final List<SAMRecord> reads, boolean isForward, final ConsensusState consensusState)
    {
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

            for(int i = 0; i < readCount; ++i)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                SAMRecord read = reads.get(i);

                locationBases[i] = 0;

                int readIndex;
                if(isForward)
                {
                    readIndex = baseIndex;

                    if(readOffsets[i] != 0 && baseIndex >= read.getReadBases().length)
                        continue;
                }
                else
                {
                    readIndex = baseIndex + readOffsets[i];

                    if(readIndex < 0)
                        continue;
                }

                locationBases[i] = reads.get(i).getReadBases()[readIndex];
                locationQuals[i] = reads.get(i).getBaseQualities()[readIndex];
                hasMismatch |= (i > 0 && locationBases[i] != locationBases[0]);
                maxQual = max(locationQuals[i], maxQual);
            }

            if(!hasMismatch)
            {
                consensusState.Bases[baseIndex] = locationBases[0];
                consensusState.BaseQualities[baseIndex] = (byte)maxQual;
            }
            else
            {
                int basePosition = consensusState.MinUnclippedPosStart + baseIndex;

                byte[] consensusBaseAndQual = determineBaseAndQual(
                        locationBases, locationQuals, reads.get(0).getContig(), basePosition);

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }
        }

        consensusState.setOutcome(ALIGNMENT_ONLY);
    }

    public byte[] determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position)
    {
        List<Byte> distinctBases = Lists.newArrayListWithCapacity(4);
        List<Integer> qualTotals = Lists.newArrayListWithCapacity(4);

        for(int i = 0; i < locationBases.length; ++i)
        {
            boolean found = false;

            for(int j = 0; j < distinctBases.size(); ++j)
            {
                if(distinctBases.get(j) == locationBases[i])
                {
                    int qualTotal = qualTotals.get(j) + locationQuals[i];
                    qualTotals.set(j, qualTotal);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                distinctBases.add(locationBases[i]);
                qualTotals.add((int)locationQuals[i]);
            }
        }

        byte maxBase = distinctBases.get(0);
        boolean maxIsRef = false;
        int maxQual = qualTotals.get(0);

        for(int i = 1; i < distinctBases.size(); ++i)
        {
            if(qualTotals.get(i) > maxQual)
            {
                maxQual = qualTotals.get(i);
                maxBase = distinctBases.get(i);
            }
            else if(qualTotals.get(i) >= maxQual && !maxIsRef)
            {
                String refBase = mRefGenome.getBaseString(chromosome, position, position);

                if(maxBase == refBase.getBytes()[0])
                {
                    maxIsRef = true;
                }
                else if(distinctBases.get(i) == refBase.getBytes()[0])
                {
                    maxQual = qualTotals.get(i);
                    maxBase = distinctBases.get(i);
                    maxIsRef = true;
                }
            }
        }

        return new byte[] { maxBase, (byte)maxQual };
    }

}
