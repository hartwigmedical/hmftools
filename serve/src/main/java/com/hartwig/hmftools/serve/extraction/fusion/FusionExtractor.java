package com.hartwig.hmftools.serve.extraction.fusion;

import static com.hartwig.hmftools.serve.extraction.fusion.FusionAnnotationConfig.EXONIC_FUSIONS_MAP;

import com.google.common.annotations.VisibleForTesting;
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
            } else {
                KnownFusionPair pair = toKnownFusionPair(event);
                if (pair == null) {
                    LOGGER.warn("Could not resolve fusion pair from '{}'", event);
                } else {
                    return validate(toKnownFusionPair(event));
                }
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
    @VisibleForTesting
    static KnownFusionPair toKnownFusionPair(@NotNull String event) {
        String[] fusionArray = event.split("-");
        String geneUp = null;
        String geneDown = null;
        if (fusionArray.length == 2) {
             geneUp = fusionArray[0];
             geneDown = fusionArray[1].split(" ")[0];
        } else if (fusionArray.length == 3) {
            // Assume one of the genes looks like "NAME-NUMBER"
            String part3 = fusionArray[2].split(" ")[0];
            if (isInteger(fusionArray[1])) {
                geneUp = fusionArray[0] + "-" + fusionArray[1];
                geneDown = part3;
            } else if (isInteger(part3)) {
                geneUp = fusionArray[0];
                geneDown = fusionArray[1] + "-" + part3;
            }
        }

        if (geneUp == null || geneDown == null) {
            return null;
        }

        return ImmutableKnownFusionPair.builder().geneUp(geneUp).geneDown(geneDown).build();
    }

    private static boolean isInteger(@NotNull String string) {
        try {
            Integer.parseInt(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
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
            if (!isIncludedSomewhereInFusionCache(pair.geneUp(), pair.geneDown())) {
                LOGGER.warn("Fusion '{}-{}' is not part of the known fusion cache", pair.geneUp(), pair.geneDown());
            }
            return pair;
        }
        return null;
    }

    private boolean isIncludedSomewhereInFusionCache(@NotNull String fiveGene, @NotNull String threeGene) {
        return knownFusionCache.hasPromiscuousFiveGene(fiveGene) || knownFusionCache.hasPromiscuousThreeGene(threeGene)
                || knownFusionCache.hasKnownFusion(fiveGene, threeGene) || knownFusionCache.hasKnownIgFusion(fiveGene, threeGene)
                || knownFusionCache.hasPromiscuousIgFusion(fiveGene);
    }
}
