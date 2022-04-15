package com.hartwig.hmftools.patientreporter.remarks;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SpecialRemarkFile {

    private static final Logger LOGGER = LogManager.getLogger(SummaryFile.class);

    private static final String SEPARATOR = "\t";

    private static final String LINE_BREAK_STRING = " <enter> ";

    private SpecialRemarkFile() {
    }

    @NotNull
    public static SpecialRemarkModel buildFromTsv(@NotNull String sampleSpecialRemarkTsv) throws IOException {
        List<String> linesSampleSpecialRemark = Files.readAllLines(new File(sampleSpecialRemarkTsv).toPath());

        Map<String, String> sampleToSpecialRemarkMap = Maps.newHashMap();

        for (String line : linesSampleSpecialRemark) {
            String[] parts = line.split(SEPARATOR);
            if (parts.length == 2) {
                String sampleId = parts[0].trim();
                String specialRemarkOfSample = parts[1].trim();
                specialRemarkOfSample = specialRemarkOfSample.replace(LINE_BREAK_STRING, "\n");
                sampleToSpecialRemarkMap.put(sampleId, specialRemarkOfSample);
            } else {
                LOGGER.warn("Suspicious line detected in sample special remark tsv: {}", line);
            }
        }

        return new SpecialRemarkModel(sampleToSpecialRemarkMap);
    }
}