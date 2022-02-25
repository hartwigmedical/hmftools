package com.hartwig.hmftools.serve.extraction.hotspot;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentMode;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HotspotExtractor.class);
    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final EventPreprocessor proteinAnnotationExtractor;
    @NotNull
    private final DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation;
    @NotNull
    private final List<DriverGene> driverGenes;

    public HotspotExtractor(@NotNull final GeneChecker geneChecker, @NotNull final ProteinResolver proteinResolver,
            @NotNull final EventPreprocessor proteinAnnotationExtractor,
            @NotNull final DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation,
            @NotNull final List<DriverGene> driverGenes) {
        this.geneChecker = geneChecker;
        this.proteinResolver = proteinResolver;
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
        this.dealWithDriverInconsistentModeAnnotation = dealWithDriverInconsistentModeAnnotation;
        this.driverGenes = driverGenes;
    }

    @Nullable
    @VisibleForTesting
    public static DriverCategory findByGene(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                return driverGene.likelihoodType();
            }
        }
        return null;
    }

    @Nullable
    public List<VariantHotspot> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.HOTSPOT && geneChecker.isValidGene(gene)) {
            DriverCategory driverCategory = findByGene(driverGenes, gene);
            if (DealWithDriverInconsistentMode.filterOnInconsistenties(dealWithDriverInconsistentModeAnnotation)) { //filter + war_only
                if (driverCategory != null) {
                    return proteinResolver.resolve(gene, transcriptId, proteinAnnotationExtractor.apply(event));
                } else {
                    if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                            DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
                        LOGGER.warn("Hotpot event on {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                        return proteinResolver.resolve(gene, transcriptId, proteinAnnotationExtractor.apply(event));
                    } else if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                            DealWithDriverInconsistentModeAnnotation.FILTER)) {
                        LOGGER.info("Hotspot event filtered -- {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                        return null;
                    }
                }
            } else { //ignore
                return proteinResolver.resolve(gene, transcriptId, proteinAnnotationExtractor.apply(event));
            }
        }
        return null;
    }
}
