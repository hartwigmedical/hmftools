package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurpleQCFile;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;

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
        List<QCValue> qcValues = Lists.newArrayList();

        List<String> qc_lines = LineReader.build().readLines(new File(path).toPath(), x -> x.contains("QCStatus"));
        assert qc_lines.size() == 1;
        String qcStatus = qc_lines.get(0).split("\t")[1];

        List<String> contamination_lines = LineReader.build().readLines(new File(path).toPath(), x -> x.contains("Contamination"));
        assert contamination_lines.size() == 1;
        String contamination = contamination_lines.get(0).split("\t")[1];

        return Lists.newArrayList(
            ImmutableQCValue.of(QCValueType.PURPLE_QC_STATUS, qcStatus),
            ImmutableQCValue.of(QCValueType.PURPLE_CONTAMINATION, contamination)
        );

    }
}
