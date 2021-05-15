package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.purity.PurpleQCFile;
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
        PurpleQC qc = PurpleQCFile.read(path);

        return Lists.newArrayList(ImmutableQCValue.builder().type(QCValueType.PURPLE_QC_STATUS).value(qc.toString()).build(),
                ImmutableQCValue.builder().type(QCValueType.PURPLE_CONTAMINATION).value(String.valueOf(qc.contamination())).build());
    }
}
