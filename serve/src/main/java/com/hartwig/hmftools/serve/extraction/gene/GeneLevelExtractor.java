package com.hartwig.hmftools.serve.extraction.gene;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentMode;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneLevelExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelExtractor.class);

    @NotNull
    private final GeneChecker exomeGeneChecker;
    @NotNull
    private final GeneChecker fusionGeneChecker;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final Set<String> activationKeyPhrases;
    @NotNull
    private final Set<String> inactivationKeyPhrases;
    private final DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation;

    public GeneLevelExtractor(@NotNull final GeneChecker exomeGeneChecker, @NotNull final GeneChecker fusionGeneChecker,
            @NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final Set<String> activationKeyPhrases, @NotNull final Set<String> inactivationKeyPhrases,
            @NotNull DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation) {
        this.exomeGeneChecker = exomeGeneChecker;
        this.fusionGeneChecker = fusionGeneChecker;
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.activationKeyPhrases = activationKeyPhrases;
        this.inactivationKeyPhrases = inactivationKeyPhrases;
        this.dealWithDriverInconsistentModeAnnotation = dealWithDriverInconsistentModeAnnotation;
    }

    @Nullable
    public GeneLevelAnnotation extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.WILD_TYPE && exomeGeneChecker.isValidGene(gene)) {
            return extractWildTypeEvents(gene, type);
        } else if (type == EventType.GENE_LEVEL && exomeGeneChecker.isValidGene(gene)) {
            return extractGeneLevelEvent(gene, event);
        } else if (type == EventType.PROMISCUOUS_FUSION && fusionGeneChecker.isValidGene(gene)) {
            return extractPromiscuousFusion(gene);
        }

        return null;
    }

    @Nullable
    GeneLevelAnnotation extractPromiscuousFusion(@NotNull String gene) {
        if (!DealWithDriverInconsistentMode.filterOnInconsistenties(dealWithDriverInconsistentModeAnnotation)) {
            if (!geneIsPresentInFusionCache(gene)) {
                if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                        DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
                    LOGGER.warn("Promiscuous fusion '{}' is not present in the known fusion cache", gene);
                    return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
                } else if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                        DealWithDriverInconsistentModeAnnotation.FILTER)) {
                    LOGGER.info("Filtered -- Promiscuous fusion '{}' is not present in the known fusion cache", gene);
                    return null;
                }
            } else {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
            }
        } else {
            if (geneIsPresentInFusionCache(gene)) {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
            } else {
                if (!geneIsPresentInFusionCache(gene)) {
                    LOGGER.warn("Filtered -- Promiscuous fusion '{}' is not present in the known fusion cache", gene);
                    return null;
                } else {
                    return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
                }
            }
        }
        return null;
    }

    @Nullable
    GeneLevelAnnotation extractWildTypeEvents(@NotNull String gene, @NotNull EventType type) {
        DriverCategory driverCategory = findByGene(driverGenes, gene);
        if (!DealWithDriverInconsistentMode.filterOnInconsistenties(dealWithDriverInconsistentModeAnnotation)) {
            if (driverCategory == null) {
                if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                        DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
                    LOGGER.warn("{} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                    return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.WILD_TYPE).build();
                } else if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                        DealWithDriverInconsistentModeAnnotation.FILTER)) {
                    LOGGER.info("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                    return null;
                } else {
                    return null;
                }
            } else {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.WILD_TYPE).build();
            }
        } else {
            if (driverCategory != null) {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.WILD_TYPE).build();
            } else {
                LOGGER.warn("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                return null;
            }
        }
    }

    @Nullable
    private static DriverCategory findByGene(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                return driverGene.likelihoodType();
            }
        }
        return null;
    }

    @Nullable
    @VisibleForTesting
    GeneLevelAnnotation extractGeneLevelEvent(@NotNull String gene, @NotNull String event) {
        int longestActivationMatchLength = -1;
        for (String keyPhrase : activationKeyPhrases) {
            if (event.contains(keyPhrase) && keyPhrase.length() > longestActivationMatchLength) {
                longestActivationMatchLength = keyPhrase.length();
            }
        }

        int longestInactivatingMatchLength = -1;
        for (String keyPhrase : inactivationKeyPhrases) {
            if (event.contains(keyPhrase) && keyPhrase.length() > longestInactivatingMatchLength) {
                longestInactivatingMatchLength = keyPhrase.length();
            }
        }

        GeneLevelEvent driverBasedEvent = determineGeneLevelEventFromDriverGenes(driverGenes, gene);
        // If we find both an activating and inactivating event, we assume the longest event is the most important.
        // This is to support cases where an activating keyphrase is a substring of inactivating keyphrase (eg "act mut" vs "inact mut".
        GeneLevelEvent result;
        if (longestActivationMatchLength > 0 || longestInactivatingMatchLength > 0) {
            result = longestActivationMatchLength >= longestInactivatingMatchLength
                    ? GeneLevelEvent.ACTIVATION
                    : GeneLevelEvent.INACTIVATION;
        } else {
            result = driverBasedEvent;
        }

        if (!DealWithDriverInconsistentMode.filterOnInconsistenties(dealWithDriverInconsistentModeAnnotation)) {
            if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                    DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
                if (driverBasedEvent == GeneLevelEvent.ANY_MUTATION) {
                    LOGGER.warn("Gene {} not present in driver catalog. {} will never be reported", gene, result);
                    return ImmutableGeneLevelAnnotation.builder().gene(gene).event(result).build();
                } else if (result != driverBasedEvent) {
                    LOGGER.warn("Mismatch in driver gene event for '{}'. Event suggests {} while driver catalog suggests {}",
                            gene,
                            result,
                            driverBasedEvent);
                    return ImmutableGeneLevelAnnotation.builder().gene(gene).event(result).build();
                }
            } else if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                    DealWithDriverInconsistentModeAnnotation.FILTER)) {
                if (driverBasedEvent == GeneLevelEvent.ANY_MUTATION) {
                    LOGGER.info("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", gene, result);
                } else if (result != driverBasedEvent) {
                    LOGGER.info("Mismatch in driver gene event for '{}'. Event suggests {} while driver catalog suggests {}",
                            gene,
                            result,
                            driverBasedEvent);
                }
                return null;
            }
        } else {
            if (driverBasedEvent != result) {
                if (driverBasedEvent == GeneLevelEvent.ANY_MUTATION) {
                    LOGGER.warn("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", gene, result);
                } else {
                    LOGGER.warn("Filtered -- Mismatch in driver gene event for '{}'. Event suggests {} while driver catalog suggests {}",
                            gene,
                            result,
                            driverBasedEvent);
                }
                return null;
            } else {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(result).build();
            }
        }
        return null;
    }

    @NotNull
    @VisibleForTesting
    static GeneLevelEvent determineGeneLevelEventFromDriverGenes(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    return GeneLevelEvent.ACTIVATION;
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    return GeneLevelEvent.INACTIVATION;
                }
            }
        }
        return GeneLevelEvent.ANY_MUTATION;
    }

    private boolean geneIsPresentInFusionCache(@NotNull String gene) {
        return knownFusionCache.hasKnownPairGene(gene) || knownFusionCache.hasPromiscuousFiveGene(gene)
                || knownFusionCache.hasPromiscuousThreeGene(gene) || knownFusionCache.hasAnyIgFusion(gene);
    }
}
