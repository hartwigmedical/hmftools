package com.hartwig.hmftools.serve.extraction.fusion;

import static com.hartwig.hmftools.serve.extraction.fusion.FusionAnnotationConfig.EXONIC_FUSIONS_MAP;
import static com.hartwig.hmftools.serve.extraction.fusion.FusionAnnotationConfig.ODDLY_NAMED_GENES_MAP;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionExtractor {

    private static final Logger LOGGER = LogManager.getLogger(FusionExtractor.class);

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public FusionExtractor(@NotNull final GeneChecker geneChecker, @NotNull final KnownFusionCache knownFusionCache) {
        this.geneChecker = geneChecker;
        this.knownFusionCache = knownFusionCache;
    }

    @Nullable
    public KnownFusionPair extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.FUSION_PAIR) {
            if (EXONIC_FUSIONS_MAP.containsKey(event)) {
                return fromConfiguredPair(EXONIC_FUSIONS_MAP.get(event), gene);
            } else if (ODDLY_NAMED_GENES_MAP.containsKey(event)) {
                return fromConfiguredPair(ODDLY_NAMED_GENES_MAP.get(event), gene);
            } else {
                String[] fusionArray = event.split("-");
                String geneUp = fusionArray[0];
                String geneDown = fusionArray[1].split(" ")[0];
                return validate(ImmutableKnownFusionPair.builder().geneUp(geneUp).geneDown(geneDown).build());
            }
        } else if (type == EventType.FUSION_PAIR_AND_EXON) {
            if (EXONIC_FUSIONS_MAP.containsKey(event)) {
                return fromConfiguredPair(EXONIC_FUSIONS_MAP.get(event), gene);
            } else {
                LOGGER.warn("Exonic fusion not configured for '{}' on '{}'", event, gene);
            }
        }
        return null;
    }

    @Nullable
    private static KnownFusionPair fromConfiguredPair(@NotNull KnownFusionPair configuredPair, @NotNull String gene) {
        KnownFusionPair pair = ImmutableKnownFusionPair.builder().from(configuredPair).build();
        if (pair.geneUp().equals(gene) || pair.geneDown().equals(gene)) {
            return pair;
        }

        LOGGER.warn("Preconfigured fusion '{}' does not match on gene level: {}", configuredPair, gene);
        return null;
    }

    @Nullable
    private KnownFusionPair validate(@NotNull KnownFusionPair pair) {
        if (geneChecker.isValidGene(pair.geneUp()) && geneChecker.isValidGene(pair.geneDown())) {
            if (!(knownFusionCache.hasPromiscuousFiveGene(pair.geneUp()) || knownFusionCache.hasPromiscuousThreeGene(pair.geneDown())
                    || knownFusionCache.hasKnownFusion(pair.geneUp(), pair.geneDown()) || knownFusionCache.hasKnownIgFusion(pair.geneUp(), pair.geneDown()) ||
                    knownFusionCache.hasPromiscuousIgFusion(pair.geneUp()))) {
                LOGGER.warn("Fusion '{}-{}' is not part of the known fusion cache", pair.geneUp(), pair.geneDown());
            }
            return pair;
        }
        return null;
    }
}
