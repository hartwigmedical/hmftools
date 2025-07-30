package com.hartwig.hmftools.cobalt;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import tech.tablesaw.api.*;

public class CobaltUtils
{
    public static Multimap<Chromosome, ReadRatio> toCommonChromosomeMap(final Table input)
    {
        Multimap<Chromosome, ReadRatio> output = ArrayListMultimap.create();
        for (String c : input.stringColumn(CobaltColumns.CHROMOSOME).unique())
        {
            if(HumanChromosome.contains(c))
            {
                Table inputFiltered = input.where(input.stringColumn(CobaltColumns.CHROMOSOME).isEqualTo(c));

                List<ReadRatio> ratios = inputFiltered.stream().map(
                        r -> ImmutableReadRatio.builder()
                            .chromosome(c)
                            .position(r.getInt(CobaltColumns.POSITION))
                            .ratio(r.getDouble(CobaltColumns.RATIO))
                            .build()).collect(Collectors.toList());

                output.putAll(HumanChromosome.fromString(c), ratios);
            }
        }
        return output;
    }

    public static CobaltRatio rowToCobaltRatio(Row row, ChromosomePositionCodec chromosomePosCodec)
    {
        return new CobaltRatio(
                chromosomePosCodec.decodeChromosome(row.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS)),
                chromosomePosCodec.decodePosition(row.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS)),
                row.getDouble(CobaltColumns.REFERENCE_READ_DEPTH),
                row.getDouble("referenceGCRatio"),
                row.getDouble(CobaltColumns.REFERENCE_GC_CONTENT),
                row.getDouble("referenceGCDiploidRatio"),
                row.getDouble(CobaltColumns.TUMOR_READ_DEPTH),
                row.getDouble("tumorGCRatio"),
                row.getDouble(CobaltColumns.TUMOR_GC_CONTENT));
    }
}
