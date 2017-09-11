package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.cosmic.CosmicModel;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfReporterData {

    @NotNull
    public abstract GeneModel geneModel();

    @NotNull
    public abstract CosmicModel cosmicModel();

    @NotNull
    public abstract DrupFilter drupFilter();

    @NotNull
    public abstract CenterModel centerModel();

    @NotNull
    public abstract String signaturePath();
}
