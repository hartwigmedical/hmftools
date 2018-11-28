package com.hartwig.hmftools.svannotation.analysis;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.VariantAnnotator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalyzer
{
    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyzer.class);

    @NotNull
    private final VariantAnnotator mGeneDataAnnotator;
    @NotNull
    private final SvDisruptionAnalyser mDisruptionAnalyser;
    @NotNull
    private final SvFusionAnalyser mFusionAnalyser;

    public StructuralVariantAnalyzer(@NotNull final VariantAnnotator geneDataAnnotator, @NotNull Set<String> disruptionGeneIDPanel,
            @NotNull final KnownFusionsModel knownFusionsModel)
    {
        mGeneDataAnnotator = geneDataAnnotator;
        mDisruptionAnalyser = new SvDisruptionAnalyser(disruptionGeneIDPanel);
        mFusionAnalyser = new SvFusionAnalyser(knownFusionsModel);
    }

    @NotNull
    public StructuralVariantAnalysis runOnVariants(@NotNull List<EnrichedStructuralVariant> variants)
    {
        List<StructuralVariantAnnotation> annotations = annotateVariants(variants);
        return runOnAnnotations(annotations);
    }

    @NotNull
    public StructuralVariantAnalysis runOnAnnotations(@NotNull List<StructuralVariantAnnotation> annotations)
    {
        List<GeneFusion> fusions = mFusionAnalyser.findFusions(annotations);
        List<GeneDisruption> disruptions = mDisruptionAnalyser.findDisruptions(annotations);

        LOGGER.debug("found {} disruptions and {} fusions", disruptions.size(), fusions.size());

        return ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
    }

    @NotNull
    public final List<StructuralVariantAnnotation> annotateVariants(@NotNull List<EnrichedStructuralVariant> variants)
    {
        return mGeneDataAnnotator.annotateVariants(variants);
    }
}
