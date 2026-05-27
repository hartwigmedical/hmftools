package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.genome.position.GenomePositions.unionOfMaps;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;

public final class PCFPositionsSupplier
{
    public static Map<Chromosome, List<PCFPosition>> createPositions(final AmberData amberData, final CobaltData cobaltData)
    {
        Map<Chromosome,List<PCFPosition>> referenceBreakPoint = cobaltData.ReferenceSegments;
        Map<Chromosome,List<PCFPosition>> tumorBreakPoints = cobaltData.TumorSegments;
        Map<Chromosome,List<PCFPosition>> tumorBAF = amberData.TumorSegments;

        PPL_LOGGER.trace("merging reference and tumor ratio break points");

        Map<Chromosome,List<PCFPosition>> chromosomePcfPositions = unionOfMaps(unionOfMaps(tumorBreakPoints, referenceBreakPoint), tumorBAF);

        /*
        Map<Chromosome,List<PCFPosition>> chromosomePcfPositions = Maps.newHashMap();

        for(Chromosome chromosome : combinedPositions.keySet())
        {
            List<PCFPosition> pcfPositions = combinedPositions.get(chromosome).stream().collect(Collectors.toList());
            chromosomePcfPositions.put(chromosome, pcfPositions);

        // ensure min and max positions do not overlap
            for(int i = 1; i < pcfPositions.size(); ++i)
            {
                PCFPosition prevPcfPosition = pcfPositions.get(i - 1);
                PCFPosition pcfPosition = pcfPositions.get(i);

                int minPosition = max(pcfPosition.minPosition(), prevPcfPosition.position() + 1);
                pcfPosition.setMinPosition(minPosition);
            }
        }
        */

        return chromosomePcfPositions;
    }
}
