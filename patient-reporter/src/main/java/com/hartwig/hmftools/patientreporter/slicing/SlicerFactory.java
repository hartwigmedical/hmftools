package com.hartwig.hmftools.patientreporter.slicing;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class SlicerFactory {

    private static final Logger LOGGER = LogManager.getLogger(SlicerFactory.class);

    private static final String BASE_RESOURCE_PATH = "bed";
    private static final String CPCT_GENE_PANEL_BED = "CPCTGenePanel.bed";
    private static final String GIAB_HIGH_CONFIDENCE_BED = "GIABHighConfidence.bed";

    private static final String FIELD_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_COLUMN = 1;
    private static final int END_COLUMN = 2;

    private SlicerFactory() {
    }

    @NotNull
    public static Slicer cpctGenePanelSlicer() throws IOException, EmptyFileException {
        return fromBedFile(Resources.getResource(BASE_RESOURCE_PATH + File.separator + CPCT_GENE_PANEL_BED).getPath());
    }

    @NotNull
    public static Slicer giabHighConfidenceSlicer() throws IOException, EmptyFileException {
        return fromBedFile(
                Resources.getResource(BASE_RESOURCE_PATH + File.separator + GIAB_HIGH_CONFIDENCE_BED).getPath());
    }

    @NotNull
    @VisibleForTesting
    static Slicer fromBedFile(@NotNull String bedFile) throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(new File(bedFile).toPath());
        final Multimap<String, GenomeRegion> regionMap = HashMultimap.create();

        String prevChromosome = null;
        GenomeRegion prevRegion = null;
        for (String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            final String chromosome = values[CHROMOSOME_COLUMN].trim();

            final long start = Long.valueOf(values[START_COLUMN].trim());
            final long end = Long.valueOf(values[END_COLUMN].trim());

            if (end < start) {
                LOGGER.warn("Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end="
                        + end);
            } else {
                final GenomeRegion region = new GenomeRegion(start, end);
                if (prevRegion != null && chromosome.equals(prevChromosome) && prevRegion.end() >= start) {
                    LOGGER.warn("BED file is not sorted, please fix! Current=" + region + ", Previous=" + prevRegion);
                } else {
                    regionMap.put(chromosome, region);
                    prevChromosome = chromosome;
                    prevRegion = region;
                }
            }
        }

        return new Slicer(regionMap);
    }
}
