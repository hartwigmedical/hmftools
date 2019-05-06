package com.hartwig.hmftools.patientreporter.summary;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SummaryFile {
    private static final Logger LOGGER = LogManager.getLogger(SummaryFile.class);
    private static final String SEPARATOR = ";";

    private SummaryFile() {
    }

    @NotNull
    public static SummaryModel buildFromCsv(@NotNull String summarySamplesCsv) throws IOException {
        List<String> linesSummarySample = LineReader.build().readLines(new File(summarySamplesCsv).toPath(), line -> line.length() > 0);

        Map<String, String> sampleToSummaryMap = Maps.newHashMap();

        for (String line : linesSummarySample) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                String sampleId = parts[0].trim();
                String summaryOfSample = parts[1].trim();
                sampleToSummaryMap.put(sampleId, summaryOfSample);
            } else {
                LOGGER.warn("Suspicious line detected in summary csv: " + line);
            }
        }

        return new SummaryModel(sampleToSummaryMap);
    }

}
