package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConstants.CANONICAL_TELOMERE_SEQUENCES;
import static com.hartwig.hmftools.telo.TeloConstants.DEFAULT_MIN_TELE_SEQ_COUNT;
import static com.hartwig.hmftools.telo.TeloConstants.DEFAULT_PARTITION_SIZE;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.commons.compress.utils.Lists;

public class TeloUtils
{
    public static boolean hasTelomericContent(final String readBases)
    {
        return hasTelomericContent(readBases, CANONICAL_TELOMERE_SEQUENCES, DEFAULT_MIN_TELE_SEQ_COUNT);
    }

    public static boolean hasTelomericContent(final String readBases, final List<String> sequences, int requiredCount)
    {
        for(String teloSeq : sequences)
        {
            int matchCount = 0;
            int seqLength = teloSeq.length();
            int matchIndex = readBases.indexOf(teloSeq);

            while(matchIndex != -1)
            {
                ++matchCount;
                matchIndex = readBases.indexOf(teloSeq, matchIndex + seqLength);

                if(matchCount >= requiredCount)
                    return true;
            }
        }

        return false;
    }

        public static List<BaseRegion> createPartitions(final TeloConfig config)
    {
        RefGenomeCoordinates refGenomeCoords = RefGenomeCoordinates.COORDS_37;

        List<BaseRegion> partitions = Lists.newArrayList();

        int partitionSize = DEFAULT_PARTITION_SIZE;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = chromosome.toString();

            if(!config.SpecificChromosomes.isEmpty() && !config.SpecificChromosomes.contains(chrStr))
                continue;

            int chromosomeLength = refGenomeCoords.lengths().get(chromosome).intValue();

            int startPos = 0;
            while(startPos < chromosomeLength)
            {
                int endPos = startPos + partitionSize - 1;

                if(endPos + partitionSize * 0.5 > chromosomeLength)
                    endPos = chromosomeLength;

                partitions.add(new BaseRegion(chrStr, startPos, endPos));

                startPos += partitionSize;
            }
        }

        return partitions;
    }
}
