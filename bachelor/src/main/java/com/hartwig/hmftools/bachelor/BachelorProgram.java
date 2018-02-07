package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;

public class BachelorProgram {
    final String name;
    final Predicate<VariantModel> vcfProcessor;
    final Predicate<GeneCopyNumber> copyNumberProcessor;
    final Predicate<HmfGenomeRegion> disruptionProcessor;

    final List<String> RequiredEffects;
    final List<String> PanelTranscripts;

    // add white and blacklist criteria for this program


    BachelorProgram(
            final String name,
            final Predicate<VariantModel> vcfProcessor,
            final Predicate<GeneCopyNumber> copyNumberProcessor,
            final Predicate<HmfGenomeRegion> disruptionProcessor,
            final List<String> requiredEffects,
            final List<String> panelTranscripts) {
        this.name = name;
        this.vcfProcessor = vcfProcessor;
        this.copyNumberProcessor = copyNumberProcessor;
        this.disruptionProcessor = disruptionProcessor;
        this.RequiredEffects = requiredEffects;
        this.PanelTranscripts = panelTranscripts;
    }
}
