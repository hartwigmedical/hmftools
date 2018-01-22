package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.genes.CosmicGeneModel;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.patientreporter.algo.MSIAnalyzer;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfReporterData {

    @NotNull
    public abstract GeneModel geneModel();

    @NotNull
    public abstract CosmicGeneModel cosmicGeneModel();

    @NotNull
    public abstract CosmicFusionModel cosmicFusionModel();

    @NotNull
    public abstract DrupFilter drupFilter();

    @NotNull
    public abstract MSIAnalyzer msiAnalyzer();
}
