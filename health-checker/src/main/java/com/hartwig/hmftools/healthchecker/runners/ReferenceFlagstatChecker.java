package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.jetbrains.annotations.NotNull;

public class ReferenceFlagstatChecker extends FileBasedHealthChecker<Flagstat> {

    public ReferenceFlagstatChecker(@NotNull String flagstatFile) throws IOException {
        super(FlagstatFile.read(flagstatFile),
                flagstat -> List.of(ImmutableQCValue.builder()
                                .type(QCValueType.REF_PROPORTION_MAPPED)
                                .value(String.valueOf(flagstat.mappedProportion()))
                                .build(),
                        ImmutableQCValue.builder()
                                .type(QCValueType.REF_PROPORTION_DUPLICATE)
                                .value(String.valueOf(flagstat.duplicateProportion()))
                                .build()));
    }
}