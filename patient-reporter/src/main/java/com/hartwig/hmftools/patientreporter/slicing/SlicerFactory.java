package com.hartwig.hmftools.patientreporter.slicing;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
    private static Slicer fromBedFile(@NotNull String bedFile) throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(new File(bedFile).toPath());
        final List<GenomeRegion> regions = Lists.newArrayList();

        GenomeRegion currentRegion = null;
        for (String line : lines) {
            currentRegion = toRegion(line, currentRegion);
            if (currentRegion != null) {
                regions.add(currentRegion);
            }
        }

        return new Slicer(regions);
    }

    @Nullable
    private static GenomeRegion toRegion(@NotNull String line, @Nullable GenomeRegion previous) {
        final String[] values = line.split(FIELD_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long start = Long.valueOf(values[START_COLUMN].trim());
        final long end = Long.valueOf(values[END_COLUMN].trim());

        if (end < start) {
            LOGGER.warn(
                    "Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end=" + end);
            return null;
        }

        final GenomeRegion region = new GenomeRegion(chromosome, start, end);
        if (previous != null && chromosome.equals(previous.chromosome()) && previous.start() > start) {
            LOGGER.warn("BED file is not sorted, please fix! Current=" + region + ", Previous=" + previous);
        }
        return region;
    }
}
