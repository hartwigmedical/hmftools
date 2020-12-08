package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.gene.ImmutableGeneLevelAnnotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneLevelExtractor {

    @NotNull
    private final GeneChecker exomeGeneChecker;
    @NotNull
    private final GeneChecker fusionGeneChecker;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final Set<String> activationKeyPhrases;
    @NotNull
    private final Set<String> inactivationKeyPhrases;

    public GeneLevelExtractor(@NotNull final GeneChecker exomeGeneChecker, @NotNull final GeneChecker fusionGeneChecker,
            @NotNull final List<DriverGene> driverGenes, @NotNull final Set<String> activationKeyPhrases,
            @NotNull final Set<String> inactivationKeyPhrases) {
        this.exomeGeneChecker = exomeGeneChecker;
        this.fusionGeneChecker = fusionGeneChecker;
        this.driverGenes = driverGenes;
        this.activationKeyPhrases = activationKeyPhrases;
        this.inactivationKeyPhrases = inactivationKeyPhrases;
    }

    @Nullable
    public GeneLevelAnnotation extract(@NotNull String gene, @NotNull EventType type, @NotNull String event) {
        if (type == EventType.GENE_LEVEL && exomeGeneChecker.isValidGene(gene)) {
            GeneLevelEvent geneLevelEvent = extractGeneLevelEvent(gene, event);
            return ImmutableGeneLevelAnnotation.builder().gene(gene).event(geneLevelEvent).build();
        } else if (type == EventType.PROMISCUOUS_FUSION && fusionGeneChecker.isValidGene(gene)) {
            return ImmutableGeneLevelAnnotation.builder().gene(gene).event(GeneLevelEvent.FUSION).build();
        }

        return null;
    }

    @VisibleForTesting
    @NotNull
    GeneLevelEvent extractGeneLevelEvent(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : activationKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return GeneLevelEvent.ACTIVATION;
            }
        }

        for (String keyPhrase : inactivationKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return GeneLevelEvent.INACTIVATION;
            }
        }

        return determineGeneLevelEventFromDriverGenes(driverGenes, gene);
    }

    @VisibleForTesting
    @NotNull
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
}
