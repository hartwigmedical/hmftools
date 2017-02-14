package com.hartwig.hmftools.patientreporter.lims;

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

public class TumorPercentages {

    private static final Logger LOGGER = LogManager.getLogger(TumorPercentages.class);

    @NotNull
    private final Map<String, Double> tumorPercentagesPerSample;

    @NotNull
    public static TumorPercentages loadFromCsv(@NotNull final String csv) throws IOException, EmptyFileException {
        final Map<String, Double> tumorPercentagesPerSample = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(csv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(",");
            assert parts.length == 2;
            tumorPercentagesPerSample.put(parts[0],
                    Double.valueOf(parts[1].replaceAll("\"", Strings.EMPTY).trim()) / 100D);
        }
        return new TumorPercentages(tumorPercentagesPerSample);
    }

    TumorPercentages(@NotNull final Map<String, Double> tumorPercentagesPerSample) {
        this.tumorPercentagesPerSample = tumorPercentagesPerSample;
    }

    public double findTumorPercentageForSample(@NotNull final String sample) {
        final Double tumorPercentage = tumorPercentagesPerSample.get(sample);
        if (tumorPercentage == null) {
            LOGGER.warn(" Could not find tumor percentage for " + sample);
            return Double.NaN;
        }
        return tumorPercentage;
    }
}
