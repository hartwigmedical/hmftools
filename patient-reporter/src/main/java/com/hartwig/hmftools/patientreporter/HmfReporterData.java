package com.hartwig.hmftools.patientreporter;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cosmic.CosmicGeneModel;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.variants.MicrosatelliteAnalyzer;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfReporterData {

    @NotNull
    public abstract GeneModel panelGeneModel();

    @NotNull
    public abstract CosmicGeneModel cosmicGeneModel();

    @NotNull
    public abstract KnownFusionsModel knownFusionsModel();

    @NotNull
    public abstract DrupFilter drupFilter();

    @NotNull
    public abstract String refGenomeFastaFileLocation();

    @NotNull
    public abstract MicrosatelliteAnalyzer microsatelliteAnalyzer();

    @NotNull
    public abstract Multimap<String, GenomeRegion> highConfidenceRegions();
}
