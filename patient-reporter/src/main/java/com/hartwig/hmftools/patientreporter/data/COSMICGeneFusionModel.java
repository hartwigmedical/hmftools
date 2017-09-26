package com.hartwig.hmftools.patientreporter.data;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class COSMICGeneFusionModel {
    @NotNull
    public abstract List<COSMICGeneFusionData> fusions();

    @NotNull
    public abstract List<COSMICPromiscuousGene> promiscuousFivePrime();

    @NotNull
    public abstract List<COSMICPromiscuousGene> promiscuousThreePrime();
}

