package com.hartwig.hmftools.patientreporter.summary;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class SummaryFile {

    private static final Logger LOGGER = LogManager.getLogger(SummaryFile.class);

    private static final String SEPARATOR = "\t";

    private static final String LINE_BREAK_STRING = " <enter> ";

    private SummaryFile() {
    }

    @NotNull
    public static SummaryModel buildFromTsv(@NotNull String sampleSummaryTsv) throws IOException {
        List<String> linesSampleSummary = LineReader.build().readLines(new File(sampleSummaryTsv).toPath(), line -> line.length() > 0);

        Map<String, String> sampleToSummaryMap = Maps.newHashMap();

        for (String line : linesSampleSummary) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                String sampleId = parts[0].trim();
                String summaryOfSample = parts[1].trim();
                summaryOfSample = summaryOfSample.replace(LINE_BREAK_STRING, "\n");
                sampleToSummaryMap.put(sampleId, summaryOfSample);
            } else {
                LOGGER.warn("Suspicious line detected in sample summary tsv: {}", line);
            }
        }

        return new SummaryModel(sampleToSummaryMap);
    }
}
