package com.hartwig.hmftools.patientreporter.slicing;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class SlicerFactory {

    private static final Logger LOGGER = LogManager.getLogger(SlicerFactory.class);

    private static final String FIELD_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_COLUMN = 1;
    private static final int END_COLUMN = 2;
    private static final int ANNOTATION_COLUMN = 3;

    private SlicerFactory() {
    }

    @NotNull
    public static Slicer fromBedFile(@NotNull String bedFile) throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(new File(bedFile).toPath());
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();

        String prevChromosome = null;
        GenomeRegion prevRegion = null;
        for (final String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            final String chromosome = values[CHROMOSOME_COLUMN].trim();
            final String annotation = values.length > ANNOTATION_COLUMN ? values[ANNOTATION_COLUMN].trim() : null;

            // KODU: BED Files are 0-based start and 1-based end, to make length simply "end - start".
            final long start = Long.valueOf(values[START_COLUMN].trim()) + 1;
            final long end = Long.valueOf(values[END_COLUMN].trim());

            if (end < start) {
                LOGGER.warn("Invalid genome region found in chromosome " + chromosome + ": start=" + start + ", end="
                        + end);
            } else {
                final GenomeRegion region = new GenomeRegion(chromosome, start, end, annotation);
                if (prevRegion != null && chromosome.equals(prevChromosome) && prevRegion.end() >= start) {
                    LOGGER.warn("BED file is not sorted, please fix! Current=" + region + ", Previous=" + prevRegion);
                } else {
                    regionMap.put(chromosome, region);
                    prevChromosome = chromosome;
                    prevRegion = region;
                }
            }
        }

        final Slicer slicer = new Slicer(regionMap);
        LOGGER.debug("Created slicer from " + bedFile + ": " + slicer.numberOfRegions() + " regions covering "
                + slicer.numberOfBases() + " bases");
        return slicer;
    }
}
