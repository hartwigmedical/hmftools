package com.hartwig.hmftools.svannotation.analysis;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.svannotation.VariantAnnotator;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalyzer
{
    @NotNull
    private final VariantAnnotator mGeneDataAnnotator;

    private SvDisruptionAnalyser mDisruptionAnalyser;
    private SvFusionAnalyser mFusionAnalyser;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyzer.class);

    public StructuralVariantAnalyzer(@NotNull final VariantAnnotator geneDataAnnotator, @NotNull Set<String> disruptionGeneIDPanel,
            @NotNull final KnownFusionsModel knownFusionsModel)
    {
        mDisruptionAnalyser = new SvDisruptionAnalyser(disruptionGeneIDPanel);
        mFusionAnalyser = new SvFusionAnalyser(knownFusionsModel);

        // TODO (KODU) See if we can replace with a filter on gene name rather than gene ensembl ID.
        mGeneDataAnnotator = geneDataAnnotator;
    }

    public final VariantAnnotator getGeneDataAnnotator() { return mGeneDataAnnotator; }

    @NotNull
    public StructuralVariantAnalysis run(final List<EnrichedStructuralVariant> variants)
    {
        LOGGER.debug("annotating variants with gene and transcript info");
        final List<StructuralVariantAnnotation> annotations = mGeneDataAnnotator.annotateVariants(variants);

        final List<StructuralVariantAnnotation> copy = Lists.newArrayList(annotations);

        LOGGER.debug("processing fusions");
        final List<GeneFusion> fusions = mFusionAnalyser.findFusions(copy);

        LOGGER.debug("processing disruptions");
        final List<GeneDisruption> disruptions = mDisruptionAnalyser.findDisruptions(copy);

        return ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
    }

    public final List<StructuralVariantAnnotation> findAnnotations(final List<EnrichedStructuralVariant> variants)
    {
        LOGGER.debug("annotating {} variants with gene and transcript info", variants.size());
        return mGeneDataAnnotator.annotateVariants(variants);
    }

    public final List<GeneFusion> findFusions(final List<StructuralVariantAnnotation> annotations)
    {
        return mFusionAnalyser.findFusions(annotations);
    }

    public final List<GeneDisruption> findDisruptions(final List<StructuralVariantAnnotation> annotations)
    {
        return mDisruptionAnalyser.findDisruptions(annotations);
    }

    // common methods
    public static boolean intronicDisruptionOnSameTranscript(final Transcript t1, final Transcript t2)
    {
        boolean sameTranscript = t1.transcriptId().equals(t2.transcriptId());
        boolean bothIntronic = t1.isIntronic() && t2.isIntronic();
        boolean sameExonUpstream = t1.exonUpstream() == t2.exonUpstream();

        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    public static boolean isUpstream(GeneAnnotation gene)
    {
        return gene.strand() * gene.orientation() > 0;
    }

    public static List<Transcript> getIntronicTranscripts(List<Transcript> transcripts)
    {
        return transcripts.stream().filter(Transcript::isIntronic).collect(Collectors.toList());
    }

}
