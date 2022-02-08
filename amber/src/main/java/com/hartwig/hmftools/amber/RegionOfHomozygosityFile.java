package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import tech.tablesaw.api.*;

import org.jetbrains.annotations.NotNull;

public final class RegionOfHomozygosityFile
{
    private static final String AMBER_EXTENSION = ".amber.homozygousregion.tsv";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull final List<RegionOfHomozygosity> regions) throws IOException
    {
        tech.tablesaw.io.csv.CsvWriteOptions writeOptions = tech.tablesaw.io.csv.CsvWriteOptions.builder(filename).separator('\t').build();
        toTable(regions).write().csv(writeOptions);
    }

    @NotNull
    static Table toTable(@NotNull final List<RegionOfHomozygosity> regions)
    {
        final Table table = tech.tablesaw.api.Table.create();

        table.addColumns(StringColumn.create("Chromosome"),
                IntColumn.create("StartPosition"),
                IntColumn.create("EndPosition"),
                IntColumn.create("SNPCount"),
                IntColumn.create("HomCount"),
                IntColumn.create("HetCount"),
                IntColumn.create("UnclearCount"),
                StringColumn.create("Filter"));

        for (RegionOfHomozygosity region : regions)
        {
            populateRow(table.appendRow(), region);
        }
        return table;
    }

    private static void populateRow(@NotNull final Row row, @NotNull final RegionOfHomozygosity region)
    {
        row.setString("Chromosome", region.getChromosome().toString());
        row.setInt("StartPosition", region.getStart());
        row.setInt("EndPosition", region.getEnd());
        row.setInt("SNPCount", region.getSnpCount());
        row.setInt("HomCount", region.getNumHomozygous());
        row.setInt("HetCount", region.getNumHeterozygous());
        row.setString("Filter", softFilter(region));
    }

    // we have several soft filters
    private static String softFilter(@NotNull final RegionOfHomozygosity region)
    {
        List<String> softFilters = new ArrayList<>();

        if (region.getLength() < AmberConstants.HOMOZYGOUS_REGION_LONG_SIZE)
            softFilters.add("minLength");

        double hetRatio = ((double)region.getNumHeterozygous()) / region.getSnpCount();

        if (hetRatio > AmberConstants.HOMOZYGOUS_REGION_MAX_HET_RATIO)
            softFilters.add("maxHetProportion");

        if (!softFilters.isEmpty())
            return String.join(";", softFilters);
        return "PASS";
    }
}
