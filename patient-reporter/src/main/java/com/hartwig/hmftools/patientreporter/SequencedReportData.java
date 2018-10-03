package com.hartwig.hmftools.patientreporter;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.enrich.CompoundEnrichment;
import com.hartwig.hmftools.patientreporter.algo.DrupActionabilityModel;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SequencedReportData {

    @NotNull
    public abstract GeneModel panelGeneModel();

    @NotNull
    public abstract CompoundEnrichment somaticVariantEnrichment();

    @NotNull
    public abstract KnownFusionsModel knownFusionsModel();

    @NotNull
    public abstract DrupActionabilityModel drupActionabilityModel();

    @NotNull
    public abstract IndexedFastaSequenceFile refGenomeFastaFile();

    @NotNull
    public abstract Multimap<String, GenomeRegion> highConfidenceRegions();
}
