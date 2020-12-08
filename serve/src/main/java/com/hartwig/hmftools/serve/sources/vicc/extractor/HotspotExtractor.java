package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotExtractor {

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final ProteinAnnotationExtractor proteinAnnotationExtractor;

    public HotspotExtractor(@NotNull final GeneChecker geneChecker, @NotNull final ProteinResolver proteinResolver,
            @NotNull final ProteinAnnotationExtractor proteinAnnotationExtractor) {
        this.geneChecker = geneChecker;
        this.proteinResolver = proteinResolver;
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
    }

    @NotNull
    public List<VariantHotspot> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.HOTSPOT && geneChecker.isValidGene(gene)) {
            return proteinResolver.resolve(gene, transcriptId, proteinAnnotationExtractor.apply(event));
        }

        return Lists.newArrayList();
    }

    @NotNull
    public ProteinAnnotationExtractor proteinAnnotationExtractor() {
        return proteinAnnotationExtractor;
    }
}
