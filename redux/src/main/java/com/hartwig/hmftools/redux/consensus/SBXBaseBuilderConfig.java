package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;

import java.util.Map;

import com.google.common.collect.Maps;

public class SBXBaseBuilderConfig extends BaseBuilderConfig
{
    public SBXBaseBuilderConfig(final RefGenome refGenome, final ConsensusStatistics consensusStats)
    {
        super(refGenome, false, false, consensusStats);
    }

    @Override
    public byte[] determineBaseAndQual(final boolean isDualStrand, final boolean[] isFirstInPair, final byte[] locationBases,
            final byte[] locationQuals, final String chromosome, final int position)
    {
        if(locationBases.length == 1)
        {
            if(locationQuals[0] == 0 && position != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, position);
                return new byte[] { refBase, 1 };
            }

            return new byte[] { locationBases[0], locationQuals[0] };
        }

        Map<Byte, int[]> baseCountsByQual = Maps.newHashMap();
        for(int i = 0; i < locationBases.length; ++i)
        {
            byte base = locationBases[i];
            if(base == NO_BASE)
            {
                continue;
            }

            int baseIdx = baseIndex(base);
            if(baseIdx < 0 || baseIdx >= DNA_BASE_BYTES.length)
            {
                continue;
            }

            byte qual = locationQuals[i];
            int[] baseCounts = baseCountsByQual.get(qual);
            if(baseCounts == null)
            {
                baseCounts = new int[] { 0, 0, 0, 0 };
                baseCountsByQual.put(qual, baseCounts);
            }

            baseCounts[baseIdx]++;
        }

        if(baseCountsByQual.isEmpty())
        {
            if(position != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, position);
                return new byte[] { refBase, 1 };
            }

            return new byte[] { NO_BASE, 0 };
        }

        if(!baseCountsByQual.containsKey((byte) SIMPLEX_QUAL) && !baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
        {
            if(position != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, position);
                return new byte[] { refBase, 1 };
            }

            return new byte[] { NO_BASE, 0 };
        }

        if(!baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
        {
            int[] simplexCounts = baseCountsByQual.get((byte) SIMPLEX_QUAL);
            int maxIdx = -1;
            int maxCount = -1;
            int totalCount = 0;
            for(int i = 0; i < simplexCounts.length; i++)
            {
                totalCount += simplexCounts[i];
                if(simplexCounts[i] > maxCount)
                {
                    maxCount = simplexCounts[i];
                    maxIdx = i;
                }
            }

            if(2 * maxCount > totalCount)
            {
                return new byte[] { DNA_BASE_BYTES[maxIdx], (byte) SIMPLEX_QUAL };
            }

            if(position != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, position);
                return new byte[] { refBase, 1 };
            }

            return new byte[] { NO_BASE, 0 };
        }

        int[] duplexCounts = baseCountsByQual.get((byte) DUPLEX_QUAL);
        int maxIdx = -1;
        boolean multipleMax = false;
        int maxCount = -1;
        int totalCount = 0;
        for(int i = 0; i < duplexCounts.length; i++)
        {
            totalCount += duplexCounts[i];
            if(duplexCounts[i] > maxCount)
            {
                maxCount = duplexCounts[i];
                maxIdx = i;
                multipleMax = false;
            }
            else if(duplexCounts[i] == maxCount)
            {
                multipleMax = true;
            }
        }

        int[] duplexErrorCounts = baseCountsByQual.get((byte) DUPLEX_ERROR_QUAL);
        if(duplexErrorCounts != null)
        {
            for(int i = 0; i < duplexErrorCounts.length; i++)
            {
                totalCount += duplexErrorCounts[i];
            }
        }

        if(2 * maxCount > totalCount)
        {
            return new byte[] { DNA_BASE_BYTES[maxIdx], (byte) DUPLEX_QUAL };
        }

        if(multipleMax)
        {
            if(position != INVALID_POSITION)
            {
                byte refBase = mRefGenome.getRefBase(chromosome, position);
                return new byte[] { refBase, 1 };
            }

            return new byte[] { NO_BASE, 0 };
        }

        return new byte[] { DNA_BASE_BYTES[maxIdx], (byte) SIMPLEX_QUAL };
    }
}
