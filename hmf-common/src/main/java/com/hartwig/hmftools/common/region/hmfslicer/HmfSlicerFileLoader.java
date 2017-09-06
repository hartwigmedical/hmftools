package com.hartwig.hmftools.common.region.hmfslicer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HmfSlicerFileLoader {
    private static final Logger LOGGER = LogManager.getLogger(BEDFileLoader.class);

    private static final String FIELD_SEPARATOR = "\t";

    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_COLUMN = 1;
    private static final int END_COLUMN = 2;
    private static final int TRANSCRIPT_ID_COLUMN = 3;
    private static final int TRANSCRIPT_VERSION_COLUMN = 4;
    private static final int GENE_COLUMN = 5;
    private static final int GENE_ID_COLUMN = 6;
    private static final int GENE_START_COLUMN = 7;
    private static final int GENE_END_COLUMN = 8;
    private static final int CHROMOSOME_BAND_COLUMN = 9;
    private static final int ENTREZ_ID_COLUMN = 10;
    private static final int EXON_ID_COLUMN = 11;
    private static final int EXON_START_COLUMN = 12;
    private static final int EXON_END_COLUMN = 13;

    private HmfSlicerFileLoader() {
    }

    @NotNull
    public static SortedSetMultimap<String, HmfGenomeRegion> fromInputStream(@NotNull final InputStream inputStream)
            throws IOException, EmptyFileException {
        return fromLines(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static SortedSetMultimap<String, HmfGenomeRegion> fromLines(@NotNull final List<String> lines)
            throws IOException, EmptyFileException {
        final SortedSetMultimap<String, HmfGenomeRegion> regionMap = TreeMultimap.create();

        String gene = "";
        ImmutableHmfGenomeRegion.Builder builder = null;

        for (final String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            final String chromosome = values[CHROMOSOME_COLUMN].trim();

            final long start = Long.valueOf(values[START_COLUMN].trim());
            final long end = Long.valueOf(values[END_COLUMN].trim());

            if (end < start) {
                LOGGER.warn("Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end=" + end);
            } else {
                if (builder == null || !gene.equals(values[GENE_COLUMN])) {
                    if (builder != null) {
                        HmfGenomeRegion region = builder.build();
                        regionMap.put(region.chromosome(), region);
                    }

                    gene = values[GENE_COLUMN];
                    builder = ImmutableHmfGenomeRegion.builder()
                            .chromosome(chromosome)
                            .start(start)
                            .end(end)
                            .transcriptID(values[TRANSCRIPT_ID_COLUMN])
                            .transcriptVersion(Integer.valueOf(values[TRANSCRIPT_VERSION_COLUMN]))
                            .chromosomeBand(values[CHROMOSOME_BAND_COLUMN])
                            .entrezId(values[ENTREZ_ID_COLUMN])
                            .gene(gene)
                            .geneID(values[GENE_ID_COLUMN])
                            .geneStart(Long.valueOf(values[GENE_START_COLUMN]))
                            .geneEnd(Long.valueOf(values[GENE_END_COLUMN]));
                }

                builder.addExome(ImmutableHmfExonRegion.builder()
                        .chromosome(chromosome)
                        .exonID(values[EXON_ID_COLUMN])
                        .start(Long.valueOf(values[EXON_START_COLUMN]))
                        .end(Long.valueOf(values[EXON_END_COLUMN]))
                        .build());
            }
        }

        if (builder != null) {
            final HmfGenomeRegion region = builder.build();
            regionMap.put(region.chromosome(), region);
        }

        return regionMap;
    }
}
