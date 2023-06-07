package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.List;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import tech.tablesaw.api.*;

public final class DiploidRatioSupplier
{
    public static Table calcDiploidRatioResults(final Table normalRatios, final List<MedianRatio> medianRatios)
    {
        Table results = normalRatios.emptyCopy();

        for (CobaltChromosome cobaltChromosome : new CobaltChromosomes(medianRatios).chromosomes())
        {
            String chr = cobaltChromosome.contig();
            if(HumanChromosome.contains(chr))
            {
                Table chrTable = normalRatios.where(normalRatios.stringColumn(CobaltColumns.CHROMOSOME).isEqualTo(chr));

                final List<Double> ratios = chrTable.doubleColumn(CobaltColumns.RATIO).asList();
                final List<Double> adjustedRatios;
                if (HumanChromosome.fromString(chr).equals(HumanChromosome._Y))
                {
                    adjustedRatios = ratios;
                }
                else
                {
                    double expectedRatio = cobaltChromosome.actualRatio();
                    adjustedRatios = new DiploidRatioNormalization(expectedRatio,
                            ROLLING_MEDIAN_MAX_DISTANCE,
                            ROLLING_MEDIAN_MIN_COVERAGE,
                            ratios).get();
                }

                chrTable.replaceColumn(DoubleColumn.create(CobaltColumns.RATIO, adjustedRatios));
                results.append(chrTable);
            }
        }

        return results;
    }
}
