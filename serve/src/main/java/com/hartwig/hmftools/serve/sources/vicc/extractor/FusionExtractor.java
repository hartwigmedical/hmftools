package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.serve.fusion.FusionAnnotationConfig.EXONIC_FUSIONS_MAP;
import static com.hartwig.hmftools.serve.fusion.FusionAnnotationConfig.ODDLY_NAMED_GENES_MAP;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(FusionExtractor.class);

    @NotNull
    private final GeneChecker geneChecker;

    public FusionExtractor(@NotNull final GeneChecker geneChecker) {
        this.geneChecker = geneChecker;
    }

    @Nullable
    public KnownFusionPair extract(@NotNull Knowledgebase source, @NotNull String gene, @NotNull EventType type, @NotNull String event) {
        ImmutableKnownFusionPair.Builder fusionBuilder = ImmutableKnownFusionPair.builder().addSources(source);

        if (type == EventType.FUSION_PAIR) {
            KnownFusionPair annotatedFusion;
            if (EXONIC_FUSIONS_MAP.containsKey(event)) {
                annotatedFusion = fusionBuilder.from(EXONIC_FUSIONS_MAP.get(event)).build();
            } else if (ODDLY_NAMED_GENES_MAP.containsKey(event)) {
                annotatedFusion = fusionBuilder.from(ODDLY_NAMED_GENES_MAP.get(event)).build();
            } else {
                String[] fusionArray = event.split("-");
                annotatedFusion = fusionBuilder.geneUp(fusionArray[0]).geneDown(fusionArray[1].split(" ")[0]).build();
            }

            if (geneChecker.isValidGene(annotatedFusion.geneUp()) && geneChecker.isValidGene(annotatedFusion.geneDown())) {
                return annotatedFusion;
            }
        } else if (type == EventType.FUSION_PAIR_AND_EXON) {
            if (EXONIC_FUSIONS_MAP.containsKey(event)) {
                KnownFusionPair fusionPair = fusionBuilder.from(EXONIC_FUSIONS_MAP.get(event)).build();
                if (fusionPair.geneUp().equals(gene) && fusionPair.geneDown().equals(gene)) {
                    return fusionPair;
                } else {
                    LOGGER.warn("Configured fusion for '{}' on '{}' does not match in terms of genes!", event, gene);
                }
            } else {
                LOGGER.warn("Exonic fusion not configured for '{}' on '{}'", event, gene);
            }
        }
        return null;
    }
}
