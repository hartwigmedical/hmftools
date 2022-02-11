package com.hartwig.hmftools.serve.extraction.codon;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonRangeByRank;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentMode;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.EnsemblFunctions;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CodonExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CodonExtractor.class);

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final MutationTypeFilterAlgo mutationTypeFilterAlgo;
    @NotNull
    private final EnsemblDataCache ensemblDataCache;
    @NotNull
    private final DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation;
    @NotNull
    private final List<DriverGene> driverGenes;

    public CodonExtractor(@NotNull final GeneChecker geneChecker, @NotNull final MutationTypeFilterAlgo mutationTypeFilterAlgo,
            @NotNull final EnsemblDataCache ensemblDataCache,
            @NotNull final DealWithDriverInconsistentModeAnnotation dealWithDriverInconsistentModeAnnotation,
            @NotNull final List<DriverGene> driverGenes) {
        this.geneChecker = geneChecker;
        this.mutationTypeFilterAlgo = mutationTypeFilterAlgo;
        this.ensemblDataCache = ensemblDataCache;
        this.dealWithDriverInconsistentModeAnnotation = dealWithDriverInconsistentModeAnnotation;
        this.driverGenes = driverGenes;
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
    public List<CodonAnnotation> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.CODON && geneChecker.isValidGene(gene)) {
            DriverCategory driverCategory = findByGene(driverGenes, gene);
            if (!DealWithDriverInconsistentMode.filterOnInconsistenties(dealWithDriverInconsistentModeAnnotation)) {
                if (driverCategory != null) {
                    if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                            DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
                        LOGGER.warn("{} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                    } else if (dealWithDriverInconsistentModeAnnotation.logging() && dealWithDriverInconsistentModeAnnotation.equals(
                            DealWithDriverInconsistentModeAnnotation.FILTER)) {
                        LOGGER.info("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                        return null;
                    }
                }
            } else {
                if (driverCategory == null) {
                    if (dealWithDriverInconsistentModeAnnotation.logging()) {
                        LOGGER.info("Filtered -- {} on {} is not included in driver catalog and won't ever be reported.", type, gene);
                        return null;
                    }
                }
            }

            HmfTranscriptRegion canonicalTranscript = EnsemblFunctions.findCanonicalTranscript(ensemblDataCache, gene);
            assert canonicalTranscript != null;

            if (transcriptId == null || transcriptId.equals(canonicalTranscript.transName())) {
                Integer codonRank = extractCodonRank(event);
                if (codonRank == null) {
                    LOGGER.warn("Could not extract codon rank from '{}'", event);
                    return null;
                }

                MutationTypeFilter mutationTypeFilter = mutationTypeFilterAlgo.determine(gene, event);
                List<CodonAnnotation> codonAnnotations =
                        determineCodonAnnotations(gene, canonicalTranscript, codonRank, mutationTypeFilter);

                if (codonAnnotations == null) {
                    LOGGER.warn("Could not resolve codon rank {} on transcript '{}' for gene '{}'",
                            codonRank,
                            canonicalTranscript.transName(),
                            gene);
                }

                return codonAnnotations;
            } else {
                LOGGER.warn("Transcript IDs not equal for provided transcript '{}' and HMF canonical transcript '{}' for {} ",
                        transcriptId,
                        canonicalTranscript.transName(),
                        event);
            }
        }

        return null;
    }

    @Nullable
    @VisibleForTesting
    static Integer extractCodonRank(@NotNull String event) {
        String codonPart;
        if (event.contains(" ")) {
            codonPart = event.split(" ")[1];
        } else {
            codonPart = event;
        }
        codonPart = codonPart.replaceAll("\\D+", "");
        if (isInteger(codonPart)) {
            return Integer.parseInt(codonPart);
        }

        return null;
    }

    @Nullable
    private static List<CodonAnnotation> determineCodonAnnotations(@NotNull String gene, @NotNull HmfTranscriptRegion canonicalTranscript,
            @Nullable Integer codonRank, @NotNull MutationTypeFilter mutationTypeFilter) {

        List<GenomeRegion> regions = Lists.newArrayList();
        if (codonRank != null) {
            regions = codonRangeByRank(canonicalTranscript, codonRank, codonRank);
        }

        if (regions != null) {
            List<CodonAnnotation> codonAnnotations = Lists.newArrayList();
            for (GenomeRegion region : regions) {
                codonAnnotations.add(ImmutableCodonAnnotation.builder()
                        .gene(gene)
                        .transcript(canonicalTranscript.transName())
                        .chromosome(region.chromosome())
                        .start(region.start())
                        .end(region.end())
                        .mutationType(mutationTypeFilter)
                        .rank(codonRank)
                        .build());
            }
            return codonAnnotations;
        }

        return null;
    }

    private static boolean isInteger(@NotNull String string) {
        try {
            Integer.parseInt(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
}
