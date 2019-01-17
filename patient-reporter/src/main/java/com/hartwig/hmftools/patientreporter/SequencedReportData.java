package com.hartwig.hmftools.patientreporter;

import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.enrich.SomaticEnrichment;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.SvAnalyzerModel;

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
    public abstract ActionabilityAnalyzer actionabilityAnalyzer();

    @NotNull
    public abstract SomaticEnrichment somaticVariantEnrichment();

    @NotNull
    public abstract KnownFusionsModel knownFusionsModel();

    @NotNull
    public abstract IndexedFastaSequenceFile refGenomeFastaFile();

    @NotNull
    public abstract Multimap<String, GenomeRegion> highConfidenceRegions();

    @NotNull
    public abstract SvAnalyzerModel svAnalyzerModel();

}
