package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FlagstatChecker implements HealthChecker {

    @NotNull
    private final String refFlagstatFile;
    @Nullable
    private final String tumFlagstatFile;

    public FlagstatChecker(@NotNull final String refFlagstatFile, @Nullable final String tumFlagstatFile) {
        this.refFlagstatFile = refFlagstatFile;
        this.tumFlagstatFile = tumFlagstatFile;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        List<QCValue> qcValues = Lists.newArrayList();

        Flagstat refFlagstat = FlagstatFile.loadFromFile(refFlagstatFile);
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_PROPORTION_MAPPED)
                .value(String.valueOf(refFlagstat.mappedProportion()))
                .build());
        qcValues.add(ImmutableQCValue.builder()
                .type(QCValueType.REF_PROPORTION_DUPLICATE)
                .value(String.valueOf(refFlagstat.duplicateProportion()))
                .build());

        if (tumFlagstatFile != null) {
            Flagstat tumFlagstat = FlagstatFile.loadFromFile(tumFlagstatFile);
            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_PROPORTION_MAPPED)
                    .value(String.valueOf(tumFlagstat.mappedProportion()))
                    .build());
            qcValues.add(ImmutableQCValue.builder()
                    .type(QCValueType.TUM_PROPORTION_DUPLICATE)
                    .value(String.valueOf(tumFlagstat.duplicateProportion()))
                    .build());
        }

        return qcValues;
    }
}
