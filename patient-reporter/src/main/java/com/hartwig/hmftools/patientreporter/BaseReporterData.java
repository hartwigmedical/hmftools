package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BaseReporterData {
    @NotNull
    public abstract CpctEcrfModel cpctEcrfModel();

    @NotNull
    public abstract LimsJsonModel limsModel();

    @NotNull
    public abstract CenterModel centerModel();

    @NotNull
    public abstract String signaturePath();
}
