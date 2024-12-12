package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.adjustBaseQual;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.MODC_BASE;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.decodeMMTag;

import org.apache.commons.lang3.NotImplementedException;

import htsjdk.samtools.SAMRecord;

public class BiomodalBaseBuilderConfig extends DefaultBaseBuilderConfig
{
    public BiomodalBaseBuilderConfig(final RefGenome refGenome, final ConsensusStatistics consensusStats)
    {
        super(refGenome, false, false, consensusStats);
    }

    @Override
    public byte[] determineBaseAndQual(final boolean isDualStrand, final boolean[] isFirstInPair, final byte[] locationBases,
            final byte[] locationQuals, final String chromosome, final int position)
    {
        int totalCQual = 0;
        int totalModCQual = 0;
        byte[] collapsedLocationBases = new byte[locationBases.length];
        for(int i = 0; i < locationBases.length; i++)
        {
            byte base = locationBases[i];
            byte qual = locationQuals[i];

            collapsedLocationBases[i] = base == MODC_BASE ? (byte) 'C' : base;

            if(base == (byte) 'C')
            {
                totalCQual += (int) qual;
            }
            else if(base == MODC_BASE)
            {
                totalModCQual += (int) qual;
            }
        }

        byte[] consensusBaseAndQual;
        consensusBaseAndQual = determineBaseAndQual(collapsedLocationBases, locationQuals, chromosome, position, (byte) 0);
        if(consensusBaseAndQual[0] == (byte) 'C' && totalModCQual > totalCQual)
            consensusBaseAndQual[0] = MODC_BASE;

        consensusBaseAndQual[1] = adjustBaseQual(consensusBaseAndQual[1]);
        return consensusBaseAndQual;
    }

    @Override
    public byte[] preprocessReadBases(final SAMRecord read)
    {
        return decodeMMTag(read, MODC_BASE);
    }
}
