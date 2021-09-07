package com.hartwig.hmftools.patientdb.clinical;

import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EcrfModels {

    @NotNull
    public abstract EcrfModel cpctModel();

    @NotNull
    public abstract EcrfModel drupModel();

    @NotNull
    public abstract WideEcrfModel wideModel();
}
