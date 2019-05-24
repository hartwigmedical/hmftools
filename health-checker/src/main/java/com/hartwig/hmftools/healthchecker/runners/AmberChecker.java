package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;

public class AmberChecker implements HealthChecker {

    @NotNull
    private final String tumorSample;
    @NotNull
    private final String amberDirectory;

    public AmberChecker(@NotNull final String tumorSample, @NotNull final String amberDirectory) {
        this.tumorSample = tumorSample;
        this.amberDirectory = amberDirectory;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        final AmberQC qcCheck = AmberQCFile.read(AmberQCFile.generateFilename(amberDirectory, tumorSample));

        final QCValue bafCheck = ImmutableQCValue.of(QCValueType.AMBER_MEAN_BAF, String.valueOf(qcCheck.meanBAF()));
        final QCValue contaminationCheck = ImmutableQCValue.of(QCValueType.AMBER_CONTAMINATION, String.valueOf(qcCheck.contamination()));

        return Lists.newArrayList(bafCheck, contaminationCheck);
    }
}
