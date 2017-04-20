package com.hartwig.hmftools.common.lims;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    private Lims() {
    }

    @NotNull
    public static LimsModel buildModelFromCsv(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        final Map<String, LimsBiopsyData> limsDataPerSample = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(",");
            if (parts.length == 4) {
                final String tumorPercentageField = parts[3].replaceAll("\"", Strings.EMPTY).trim();
                // KODU: Only load tumor samples into LIMS datamodel for now.
                if (!tumorPercentageField.equals(Strings.EMPTY)) {
                    final String sample = parts[0].trim();
                    final String samplingDate = !parts[1].trim().equals(Strings.EMPTY) ? parts[1].trim() : null;
                    final String arrivalDate = parts[2].trim();
                    final double tumorPercentage = Double.valueOf(tumorPercentageField) / 100D;
                    limsDataPerSample.put(sample, new LimsBiopsyData(samplingDate, arrivalDate, tumorPercentage));
                } else {
                    LOGGER.warn("Invalid row found in lims csv: " + line);
                }
            }
        }
        return new LimsModel(limsDataPerSample);
    }

    @NotNull
    public static LimsModel buildEmptyModel() {
        return new LimsModel(Maps.newHashMap());
    }
}
