package com.hartwig.hmftools.serve.extraction.fusion;

import static com.hartwig.hmftools.serve.extraction.fusion.FusionAnnotationConfig.EXONIC_FUSIONS_MAP;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentMode;
import com.hartwig.hmftools.serve.extraction.catalog.DriverInconsistencyMode;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.apache.commons.compress.utils.Lists;
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
    @NotNull
    private final Set<String> exonicDelDupFusionKeyPhrases;
    private final DriverInconsistencyMode driverInconsistencyMode;

    public FusionExtractor(@NotNull final GeneChecker geneChecker, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final Set<String> exonicDelDupFusionKeyPhrases,
            @NotNull final DriverInconsistencyMode driverInconsistencyMode) {
        this.geneChecker = geneChecker;
        this.knownFusionCache = knownFusionCache;
        this.exonicDelDupFusionKeyPhrases = exonicDelDupFusionKeyPhrases;
        this.driverInconsistencyMode = driverInconsistencyMode;
    }

    @Nullable
    public KnownFusionPair extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.FUSION_PAIR) {
            if (EXONIC_FUSIONS_MAP.containsKey(event)) {
                return fromConfiguredPair(EXONIC_FUSIONS_MAP.get(event), gene);
            } else if (hasExonicDelDupKeyPhrase(event)) {
                return validate(fromExonicDelDup(gene, event));
            } else {
                return validate(fromStandardFusionPairEvent(event));
            }
        } else if (type == EventType.FUSION_PAIR_AND_EXON) {
            return validate(fromExonicDelDup(gene, event));
        } else {
            return null;
        }
    }

    private boolean hasExonicDelDupKeyPhrase(@NotNull String event) {
        for (String keyPhrase : exonicDelDupFusionKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }

    @Nullable
    private static KnownFusionPair fromExonicDelDup(@NotNull String gene, @NotNull String event) {
        ExonicDelDupType exonicDelDupType = FusionAnnotationConfig.DEL_DUP_TYPE_PER_GENE.get(gene);
        if (exonicDelDupType == null) {
            LOGGER.warn("No exonic del dup type configured for gene '{}'", gene);
            return null;
        }

        Integer exonRank = extractExonRank(event);
        if (exonRank == null) {
            return null;
        }

        int exonUp;
        int exonDown;
        switch (exonicDelDupType) {
            case FULL_EXONIC_DELETION: {
                exonUp = exonRank - 1;
                exonDown = exonRank + 1;
                break;
            }
            case PARTIAL_EXONIC_DELETION: {
                exonUp = exonRank;
                exonDown = exonRank;
                break;
            }
            default: {
                throw new IllegalStateException("Unrecognized del dup type: " + exonicDelDupType);
            }
        }

        return ImmutableKnownFusionPair.builder()
                .geneUp(gene)
                .minExonUp(exonUp)
                .maxExonUp(exonUp)
                .geneDown(gene)
                .minExonDown(exonDown)
                .maxExonDown(exonDown)
                .build();
    }

    @Nullable
    private static Integer extractExonRank(@NotNull String event) {
        List<Integer> exons = Lists.newArrayList();
        String[] words = event.split(" ");
        for (String word : words) {
            if (isInteger(word)) {
                exons.add(Integer.valueOf(word));
            }
        }
        if (exons.size() > 1) {
            LOGGER.warn("Multiple exon ranks extracted from '{}' while expecting 1", event);
            return null;
        } else if (exons.isEmpty()) {
            LOGGER.warn("No exon rank could be resolved from '{}'", event);
            return null;
        }

        return exons.get(0);
    }

    @Nullable
    private KnownFusionPair fromStandardFusionPairEvent(@NotNull String event) {
        String[] fusionArray = event.split("-");
        String geneUp = null;
        String geneDown = null;
        if (fusionArray.length == 2) {
            geneUp = fusionArray[0];
            geneDown = fusionArray[1].split(" ")[0];
        } else if (fusionArray.length == 3) {
            String geneUpScenario1 = fusionArray[0] + "-" + fusionArray[1];
            String geneDownScenario1 = fusionArray[2].split(" ")[0];
            String geneUpScenario2 = fusionArray[0];
            String geneDownScenario2 = fusionArray[1] + "-" + fusionArray[2].split(" ")[0];

            if (geneChecker.geneExistsInAllValidGenes(geneUpScenario1) && geneChecker.geneExistsInAllValidGenes(geneDownScenario1)) {
                geneUp = geneUpScenario1;
                geneDown = geneDownScenario1;
            } else if (geneChecker.geneExistsInAllValidGenes(geneUpScenario2) && geneChecker.geneExistsInAllValidGenes(geneDownScenario2)) {
                geneUp = geneUpScenario2;
                geneDown = geneDownScenario2;
            }
        } else if (fusionArray.length == 4) {
            geneUp = fusionArray[0] + "-" + fusionArray[1];
            geneDown = fusionArray[2] + "-" + fusionArray[3].split(" ")[0];
        }

        if (geneUp == null || geneDown == null) {
            LOGGER.warn("Could not resolve fusion pair from '{}'", event);
            return null;
        }

        return ImmutableKnownFusionPair.builder().geneUp(removeAllSpaces(geneUp)).geneDown(removeAllSpaces(geneDown)).build();
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
    private KnownFusionPair fromConfiguredPair(@NotNull KnownFusionPair configuredPair, @NotNull String gene) {
        KnownFusionPair pair = ImmutableKnownFusionPair.builder().from(configuredPair).build();
        KnownFusionPair pairValidated = validate(pair);
        if (pairValidated != null) {
            if (!pairValidated.geneUp().equals(gene) || !pairValidated.geneDown().equals(gene)) {
                LOGGER.warn("Preconfigured fusion '{}' does not match on gene level: {}", configuredPair, gene);
                return null;
            }
        }

        return pair;
    }

    @Nullable
    private KnownFusionPair validate(@Nullable KnownFusionPair pair) {
        if (pair == null) {
            return null;
        }

        if (geneChecker.isValidGene(pair.geneUp()) && geneChecker.isValidGene(pair.geneDown())) {
            if (DealWithDriverInconsistentMode.filterOnInconsistencies(driverInconsistencyMode)) {
                if (!isIncludedSomewhereInFusionCache(pair.geneUp(), pair.geneDown())) {
                    if (driverInconsistencyMode.logging() && driverInconsistencyMode.equals(
                            DriverInconsistencyMode.WARN_ONLY)) {
                        LOGGER.warn("Fusion event on fusion '{}-{}' is not part of the known fusion cache", pair.geneUp(), pair.geneDown());
                        return pair;
                    } else if (driverInconsistencyMode.logging() && driverInconsistencyMode.equals(
                            DriverInconsistencyMode.FILTER)) {
                        LOGGER.info("Fusion evnet filtered -- Fusion '{}-{}' is not part of the known fusion cache", pair.geneUp(), pair.geneDown());
                        return null;
                    }
                } else {
                    return pair;
                }
            } else {
                return pair;

            }
        }
        return null;
    }

    private boolean isIncludedSomewhereInFusionCache(@NotNull String fiveGene, @NotNull String threeGene) {
        return knownFusionCache.hasExonDelDup(fiveGene) || knownFusionCache.hasPromiscuousFiveGene(fiveGene)
                || knownFusionCache.hasPromiscuousThreeGene(threeGene) || knownFusionCache.hasKnownFusion(fiveGene, threeGene)
                || knownFusionCache.hasKnownIgFusion(fiveGene, threeGene) || knownFusionCache.hasPromiscuousIgFusion(fiveGene);
    }

    @NotNull
    private static String removeAllSpaces(@NotNull String value) {
        return value.replaceAll("\\s+", "");
    }
}