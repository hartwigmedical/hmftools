package com.hartwig.hmftools.common.region.hmfslicer;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public abstract class HmfSlicerFileLoader {
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

    private HmfSlicerFileLoader() {
    }

    @NotNull
    public static SortedSetMultimap<String, HmfGenomeRegion> fromHmfGenePanelFile(@NotNull String genePanelFile)
            throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(new File(genePanelFile).toPath());
        final SortedSetMultimap<String, HmfGenomeRegion> regionMap = TreeMultimap.create();

        for (final String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            final String chromosome = values[CHROMOSOME_COLUMN].trim();

            // KODU: Positions from BED files are 0-based start and 1-based end, to make length simply "end - start".
            final long start = Long.valueOf(values[START_COLUMN].trim()) + 1;
            final long end = Long.valueOf(values[END_COLUMN].trim());

            if (end < start) {
                LOGGER.warn("Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end="
                        + end);
            } else {
                final String transcriptId = values[TRANSCRIPT_ID_COLUMN];
                final int transcriptVersion = Integer.valueOf(values[TRANSCRIPT_VERSION_COLUMN]);
                final String gene = values[GENE_COLUMN];
                final String chromosomeBand = values[CHROMOSOME_BAND_COLUMN];
                final String entrezId = values[ENTREZ_ID_COLUMN];
                final String geneId = values[GENE_ID_COLUMN];
                final long geneStart = Long.valueOf(values[GENE_START_COLUMN]);
                final long geneEnd = Long.valueOf(values[GENE_END_COLUMN]);

                final HmfGenomeRegion region = new ImmutableHmfGenomeRegion(chromosome, start, end, transcriptId,
                        transcriptVersion, gene, geneId, geneStart, geneEnd, chromosomeBand, entrezId);

                regionMap.put(chromosome, region);
            }
        }

        return regionMap;
    }
}
