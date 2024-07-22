package com.hartwig.hmftools.purple.segment;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.position.GenomePositions.union;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;

public final class PCFPositionsSupplier
{
    public static Map<Chromosome, List<PCFPosition>> createPositions(final AmberData amberData, final CobaltData cobaltData)
    {
        final Multimap<Chromosome,PCFPosition> referenceBreakPoint = cobaltData.ReferenceSegments;
        final Multimap<Chromosome,PCFPosition> tumorBreakPoints = cobaltData.TumorSegments;
        final Multimap<Chromosome,PCFPosition> tumorBAF = amberData.TumorSegments;

        PPL_LOGGER.info("merging reference and tumor ratio break points");

        Multimap<Chromosome,PCFPosition> combinedPositions = union(union(tumorBreakPoints, referenceBreakPoint), tumorBAF);

        Map<Chromosome,List<PCFPosition>> chromosomePcfPositions = Maps.newHashMap();

        for(Chromosome chromosome : combinedPositions.keySet())
        {
            List<PCFPosition> pcfPositions = combinedPositions.get(chromosome).stream().collect(Collectors.toList());
            chromosomePcfPositions.put(chromosome, pcfPositions);

        /*
        // ensure min and max positions do not overlap
            for(int i = 1; i < pcfPositions.size(); ++i)
            {
                PCFPosition prevPcfPosition = pcfPositions.get(i - 1);
                PCFPosition pcfPosition = pcfPositions.get(i);

                int minPosition = max(pcfPosition.minPosition(), prevPcfPosition.position() + 1);
                pcfPosition.setMinPosition(minPosition);
            }
        */
        }

        return chromosomePcfPositions;
    }
}
