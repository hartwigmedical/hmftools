package com.hartwig.hmftools.serve.extraction.gene;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.util.DriverInconsistencyMode;
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
    @NotNull
    private final Set<String> genericKeyPhrases;
    private final DriverInconsistencyMode driverInconsistencyMode;

    public GeneLevelExtractor(@NotNull final GeneChecker exomeGeneChecker, @NotNull final GeneChecker fusionGeneChecker,
            @NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final Set<String> activationKeyPhrases, @NotNull final Set<String> inactivationKeyPhrases,
            @NotNull final Set<String> genericKeyPhrases, @NotNull DriverInconsistencyMode driverInconsistencyMode) {
        this.exomeGeneChecker = exomeGeneChecker;
        this.fusionGeneChecker = fusionGeneChecker;
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.activationKeyPhrases = activationKeyPhrases;
        this.inactivationKeyPhrases = inactivationKeyPhrases;
        this.driverInconsistencyMode = driverInconsistencyMode;
        this.genericKeyPhrases = genericKeyPhrases;
    }

    @Nullable
    public GeneLevelAnnotation extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.WILD_TYPE && exomeGeneChecker.isValidGene(gene)) {
            return extractWildTypeEvent(gene, type);
        } else if (type == EventType.GENE_LEVEL && exomeGeneChecker.isValidGene(gene)) {
            return extractGeneLevelEvent(gene, event);
        } else if (type == EventType.PROMISCUOUS_FUSION && fusionGeneChecker.isValidGene(gene)) {
            return extractPromiscuousFusion(gene);
        }

        return null;
    }

    @Nullable
    GeneLevelAnnotation extractPromiscuousFusion(@NotNull String gene) {
        if (driverInconsistencyMode.isActive() && !geneIsPresentInFusionCache(gene)) {
            if (driverInconsistencyMode == DriverInconsistencyMode.WARN_ONLY) {
                LOGGER.warn("Promiscuous fusion '{}' is not present in the known fusion cache", gene);
            } else if (driverInconsistencyMode == DriverInconsistencyMode.FILTER) {
                LOGGER.info("Promiscuous fusion filtered -- Promiscuous fusion '{}' is not present in the known fusion cache", gene);
                return null;
            }
        }

        return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
    }

    @Nullable
    GeneLevelAnnotation extractWildTypeEvent(@NotNull String gene, @NotNull EventType type) {
        DriverCategory driverCategory = findByGene(driverGenes, gene);

        if (driverCategory == null && driverInconsistencyMode.isActive()) {
            if (driverInconsistencyMode == DriverInconsistencyMode.WARN_ONLY) {
                LOGGER.warn("Wildtype event {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
            } else if (driverInconsistencyMode == DriverInconsistencyMode.FILTER) {
                LOGGER.info("Wildtype event filtered -- {} on {} is not included in driver catalog and won't ever be reported.",
                        type,
                        gene);
                return null;
            }
        }

        return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.WILD_TYPE).build();
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

    static boolean geneInDriverGenes(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                return true;
            }
        }
        return false;
    }

    @Nullable
    @VisibleForTesting
    GeneLevelAnnotation extractGeneLevelEvent(@NotNull String gene, @NotNull String event) {
        GeneLevelEvent result = GeneLevelEvent.ANY_MUTATION;
        for (String keyPhrase : activationKeyPhrases) {
            if (event.contains(keyPhrase)) {
                result = GeneLevelEvent.ACTIVATION;
            }
        }

        for (String keyPhrase : inactivationKeyPhrases) {
            if (event.contains(keyPhrase)) {
                result = GeneLevelEvent.INACTIVATION;
            }
        }

        for (String keyPhrase : genericKeyPhrases) {
            if (event.contains(keyPhrase)) {
                result = GeneLevelEvent.ANY_MUTATION;
            }
        }

        GeneLevelEvent driverBasedEvent = determineGeneLevelEventFromDriverGenes(driverGenes, gene);

        if (driverInconsistencyMode.isActive()) {
            if (driverInconsistencyMode == DriverInconsistencyMode.WARN_ONLY) {
                if (!geneInDriverGenes(driverGenes, gene)) {
                    LOGGER.warn("Gene level event on gene {} not present in driver catalog. {} will never be reported", gene, result);
                } else if (geneInDriverGenes(driverGenes, gene) && result != GeneLevelEvent.ANY_MUTATION && result != driverBasedEvent) {
                    LOGGER.warn(
                            "Gene level event mismatch in driver gene event for '{}'. Event suggests {} while driver catalog suggests {}",
                            gene,
                            result,
                            driverBasedEvent);
                }
            } else if (driverInconsistencyMode == DriverInconsistencyMode.FILTER) {
                if (!geneInDriverGenes(driverGenes, gene)) {
                    LOGGER.info("Gene level event filtered -- {} on {} is not included in driver catalog and won't ever be reported.",
                            result,
                            gene);
                    return null;
                } else if (geneInDriverGenes(driverGenes, gene) && result != GeneLevelEvent.ANY_MUTATION && result != driverBasedEvent) {
                    LOGGER.info(
                            "Gene level event filtered -- Mismatch in driver gene event for '{}'. "
                                    + "Event suggests {} while driver catalog suggests {}",
                            gene,
                            result,
                            driverBasedEvent);
                    return null;
                }
            }
        }

        return ImmutableGeneLevelAnnotation.builder().gene(gene).event(result).build();
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