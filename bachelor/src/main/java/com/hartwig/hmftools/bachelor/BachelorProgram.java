package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

public class BachelorProgram {
    private final String name;
    private final Predicate<VariantModel> vcfProcessor;
    private final Predicate<GeneCopyNumber> copyNumberProcessor;
    private final Predicate<HmfGenomeRegion> disruptionProcessor;

    private final List<String> RequiredEffects;
    private final List<String> PanelTranscripts;

    // add white and blacklist criteria for this program

    BachelorProgram(final String name, final Predicate<VariantModel> vcfProcessor, final Predicate<GeneCopyNumber> copyNumberProcessor,
            final Predicate<HmfGenomeRegion> disruptionProcessor, final List<String> requiredEffects, final List<String> panelTranscripts) {
        this.name = name;
        this.vcfProcessor = vcfProcessor;
        this.copyNumberProcessor = copyNumberProcessor;
        this.disruptionProcessor = disruptionProcessor;
        this.RequiredEffects = requiredEffects;
        this.PanelTranscripts = panelTranscripts;
    }

    public String name() {
        return name;
    }

    public Predicate<VariantModel> vcfProcessor() {
        return vcfProcessor;
    }

    public Predicate<GeneCopyNumber> copyNumberProcessor() {
        return copyNumberProcessor;
    }

    public Predicate<HmfGenomeRegion> disruptionProcessor() {
        return disruptionProcessor;
    }

    public List<String> requiredEffects() {
        return RequiredEffects;
    }

    public List<String> panelTranscripts() {
        return PanelTranscripts;
    }
}
