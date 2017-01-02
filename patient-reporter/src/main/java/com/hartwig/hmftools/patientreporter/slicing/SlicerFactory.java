package com.hartwig.hmftools.patientreporter.slicing;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.jetbrains.annotations.NotNull;

public final class SlicerFactory {

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
        List<String> lines = FileReader.build().readLines(new File(bedFile).toPath());
        List<GenomeRegion> regions = lines.stream().map(SlicerFactory::toRegion).collect(Collectors.toList());
        return new Slicer(regions);
    }

    @NotNull
    private static GenomeRegion toRegion(@NotNull String line) {
        String[] values = line.split(FIELD_SEPARATOR);

        return new GenomeRegion(values[CHROMOSOME_COLUMN], Long.valueOf(values[START_COLUMN].trim()),
                Long.valueOf(values[END_COLUMN].trim()));
    }
}
