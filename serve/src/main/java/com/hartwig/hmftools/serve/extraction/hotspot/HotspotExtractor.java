package com.hartwig.hmftools.serve.extraction.hotspot;

import java.util.List;

import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotExtractor {

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final EventPreprocessor proteinAnnotationExtractor;

    public HotspotExtractor(@NotNull final GeneChecker geneChecker, @NotNull final ProteinResolver proteinResolver,
            @NotNull final EventPreprocessor proteinAnnotationExtractor) {
        this.geneChecker = geneChecker;
        this.proteinResolver = proteinResolver;
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
    }

    @Nullable
    public List<VariantHotspot> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.HOTSPOT && geneChecker.isValidGene(gene)) {
            return proteinResolver.resolve(gene, transcriptId, proteinAnnotationExtractor.apply(event));
        }

        return null;
    }
}
