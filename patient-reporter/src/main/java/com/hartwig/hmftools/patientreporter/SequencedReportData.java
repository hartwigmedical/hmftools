package com.hartwig.hmftools.patientreporter;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

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
    public abstract IndexedFastaSequenceFile refGenomeFastaFile();

    @NotNull
    public abstract Multimap<String, GenomeRegion> highConfidenceRegions();

    @NotNull
    public abstract GermlineReportingModel germlineReportingModel();

    @NotNull
    public abstract SummaryModel summaryModel();

}
