package com.hartwig.hmftools.serve.extraction.codon;

import static com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils.codonRangeByIndex;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventType;
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
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public CodonExtractor(@NotNull final GeneChecker geneChecker, @NotNull final MutationTypeFilterAlgo mutationTypeFilterAlgo,
            @NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.geneChecker = geneChecker;
        this.mutationTypeFilterAlgo = mutationTypeFilterAlgo;
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @Nullable
    public List<CodonAnnotation> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.CODON && geneChecker.isValidGene(gene)) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(gene);
            assert canonicalTranscript != null;

            if (transcriptId == null || transcriptId.equals(canonicalTranscript.transName())) {
                Integer codonIndex = extractCodonIndex(event);
                if (codonIndex == null) {
                    LOGGER.warn("Could not extract codon index from '{}'", event);
                    return null;
                }

                MutationTypeFilter mutationTypeFilter = mutationTypeFilterAlgo.determine(gene, event);
                List<CodonAnnotation> codonAnnotations = determineCodonAnnotations(gene,
                        canonicalTranscript,
                        codonIndex,
                        mutationTypeFilter);

                if (codonAnnotations == null) {
                    LOGGER.warn("Could not resolve codon index {} on transcript '{}' for gene '{}'",
                            codonIndex,
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
    static Integer extractCodonIndex(@NotNull String event) {
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
            int codonIndex, @NotNull MutationTypeFilter mutationTypeFilter) {
        List<GenomeRegion> regions = codonRangeByIndex(canonicalTranscript, codonIndex, codonIndex);

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
                        .codonIndex(codonIndex)
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
