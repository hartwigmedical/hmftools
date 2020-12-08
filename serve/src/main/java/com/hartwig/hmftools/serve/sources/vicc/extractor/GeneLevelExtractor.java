package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.gene.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.vicc.annotation.ViccClassificationConfig;

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

    public GeneLevelExtractor(@NotNull final GeneChecker exomeGeneChecker, @NotNull final GeneChecker fusionGeneChecker,
            @NotNull final List<DriverGene> driverGenes) {
        this.exomeGeneChecker = exomeGeneChecker;
        this.fusionGeneChecker = fusionGeneChecker;
        this.driverGenes = driverGenes;
    }

    @Nullable
    public GeneLevelAnnotation extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.GENE_LEVEL && exomeGeneChecker.isValidGene(gene)) {
            GeneLevelEvent geneLevelEvent = extractGeneLevelEvent(gene, driverGenes, event);
            if (event != null) {
                return ImmutableGeneLevelAnnotation.builder().gene(gene).event(geneLevelEvent).build();
            }
        } else if (type == EventType.PROMISCUOUS_FUSION && fusionGeneChecker.isValidGene(gene)) {
            return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
        }

        return null;
    }

    @Nullable
    @VisibleForTesting
    static GeneLevelEvent extractGeneLevelEvent(@NotNull String gene, @NotNull List<DriverGene> driverGenes, @NotNull String event) {
        String trimmedEvent = event.trim();

        if (trimmedEvent.contains(" ")) {
            String firstWord = trimmedEvent.split(" ")[0];
            if (firstWord.equals(gene)) {
                trimmedEvent = trimmedEvent.split(" ", 2)[1].trim();
            }
        }

        if (ViccClassificationConfig.INACTIVATING_GENE_LEVEL_KEY_PHRASES.contains(trimmedEvent)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (ViccClassificationConfig.ACTIVATING_GENE_LEVEL_KEY_PHRASES.contains(trimmedEvent)) {
            return GeneLevelEvent.ACTIVATION;
        } else if (ViccClassificationConfig.GENERIC_GENE_LEVEL_KEY_PHRASES.contains(trimmedEvent) || gene.equals(trimmedEvent.replaceAll(
                "\\s+",
                ""))) {
            return determineGeneLevelEventFromDriverGenes(gene, driverGenes);
        }

        LOGGER.warn("Could not determine gene level event for '{}' on '{}'", event, gene);
        return null;
    }

    @VisibleForTesting
    @NotNull
    static GeneLevelEvent determineGeneLevelEventFromDriverGenes(@NotNull String gene, @NotNull List<DriverGene> driverGenes) {
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
}
