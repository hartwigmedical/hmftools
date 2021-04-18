package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurpleQCFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleChecker implements HealthChecker {

    @NotNull
    private final String tumorSample;
    @NotNull
    private final String purpleDirectory;

    public PurpleChecker(@NotNull final String tumorSample, @NotNull final String purpleDirectory) {
        this.tumorSample = tumorSample;
        this.purpleDirectory = purpleDirectory;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        String path = PurpleQCFile.generateFilename(purpleDirectory, tumorSample);

        List<String> purpleQcLines = Files.readAllLines(new File(path).toPath());
        String qcStatus = valueBySubstring(purpleQcLines, "QCStatus");
        String contamination = valueBySubstring(purpleQcLines, "Contamination");

        if (qcStatus == null || contamination == null) {
            throw new IOException("Unable to parse purple QC file correctly");
        }

        return Lists.newArrayList(ImmutableQCValue.builder().type(QCValueType.PURPLE_QC_STATUS).value(qcStatus).build(),
                ImmutableQCValue.builder().type(QCValueType.PURPLE_CONTAMINATION).value(contamination).build());
    }

    @Nullable
    private static String valueBySubstring(@NotNull List<String> lines, @NotNull String subString) {
        List<String> matchLines = Lists.newArrayList();
        for (String line : lines) {
            if (line.contains(subString)) {
                matchLines.add(line);
            }
        }
        if (matchLines.size() == 1) {
            return matchLines.get(0).split("\t")[1];
        }
        return null;
    }
}
