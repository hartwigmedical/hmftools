package com.hartwig.hmftools.cobalt;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import tech.tablesaw.api.*;
import tech.tablesaw.columns.Column;

public class CobaltUtils
{
    // helper function to set the correct column name
    public static Table replaceColumn(Table table, String columnName, Column<?> newColumn)
    {
        return table.replaceColumn(columnName, newColumn.setName(columnName));
    }

    public static Multimap<com.hartwig.hmftools.common.genome.chromosome.Chromosome, ReadRatio> toCommonChromosomeMap(
            Table input)
    {
        Multimap<com.hartwig.hmftools.common.genome.chromosome.Chromosome, ReadRatio> output = ArrayListMultimap.create();
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
        return ImmutableCobaltRatio.builder()
                .chromosome(chromosomePosCodec.decodeChromosome(row.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS)))
                .position(chromosomePosCodec.decodePosition(row.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS)))
                .referenceReadDepth(row.getDouble(CobaltColumns.REFERENCE_READ_DEPTH))
                .tumorReadDepth(row.getDouble(CobaltColumns.TUMOR_READ_DEPTH))
                .referenceGCRatio(row.getDouble("referenceGCRatio"))
                .tumorGCRatio(row.getDouble("tumorGCRatio"))
                .referenceGCDiploidRatio(row.getDouble("referenceGCDiploidRatio"))
                .referenceGcContent(row.getDouble(CobaltColumns.REFERENCE_GC_CONTENT))
                .tumorGcContent(row.getDouble(CobaltColumns.TUMOR_GC_CONTENT))
                .build();
    }

    public static Table createReadCountTable()
    {
        return Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.READ_DEPTH));
    }

    public static Table createRatioTable()
    {
        return Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.RATIO));
    }
}
