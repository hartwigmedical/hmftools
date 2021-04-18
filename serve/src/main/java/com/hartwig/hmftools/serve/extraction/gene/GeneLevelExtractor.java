package com.hartwig.hmftools.serve.extraction.gene;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventType;
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
    private final boolean reportOnDriverInconsistencies;

    public GeneLevelExtractor(@NotNull final GeneChecker exomeGeneChecker, @NotNull final GeneChecker fusionGeneChecker,
            @NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final Set<String> activationKeyPhrases, @NotNull final Set<String> inactivationKeyPhrases,
            final boolean reportOnDriverInconsistencies) {
        this.exomeGeneChecker = exomeGeneChecker;
        this.fusionGeneChecker = fusionGeneChecker;
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.activationKeyPhrases = activationKeyPhrases;
        this.inactivationKeyPhrases = inactivationKeyPhrases;
        this.reportOnDriverInconsistencies = reportOnDriverInconsistencies;
    }

    @Nullable
    public GeneLevelAnnotation extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.GENE_LEVEL && exomeGeneChecker.isValidGene(gene)) {
            GeneLevelEvent geneLevelEvent = extractGeneLevelEvent(gene, event);
            return ImmutableGeneLevelAnnotation.builder().gene(gene).event(geneLevelEvent).build();
        } else if (type == EventType.PROMISCUOUS_FUSION && fusionGeneChecker.isValidGene(gene)) {
            if (reportOnDriverInconsistencies && !geneIsPresentInFusionCache(gene)) {
                LOGGER.warn("Promiscuous fusion '{}' is not present in the known fusion cache", gene);
            }
            return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
        }

        return null;
    }

    @NotNull
    @VisibleForTesting
    GeneLevelEvent extractGeneLevelEvent(@NotNull String gene, @NotNull String event) {
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
        if (longestActivationMatchLength > 0 || longestInactivatingMatchLength > 0) {
            GeneLevelEvent geneLevelEvent = longestActivationMatchLength >= longestInactivatingMatchLength
                    ? GeneLevelEvent.ACTIVATION
                    : GeneLevelEvent.INACTIVATION;
            if (reportOnDriverInconsistencies && geneLevelEvent != driverBasedEvent && driverBasedEvent != GeneLevelEvent.ANY_MUTATION) {
                LOGGER.warn("Mismatch in driver gene event for '{}'. Event suggests {} while driver catalog suggests {}",
                        gene,
                        geneLevelEvent,
                        driverBasedEvent);
            }
            return geneLevelEvent;
        } else {
            return driverBasedEvent;
        }
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
