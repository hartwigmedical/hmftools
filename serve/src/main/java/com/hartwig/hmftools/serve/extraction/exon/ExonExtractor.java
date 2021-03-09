package com.hartwig.hmftools.serve.extraction.exon;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilterAlgo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ExonExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ExonExtractor.class);

    private static final int SPLICE_SIZE = 10;

    private static final Set<EventType> EXON_EVENTS = Sets.newHashSet(EventType.EXON, EventType.FUSION_PAIR_AND_EXON);

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final MutationTypeFilterAlgo mutationTypeFilterAlgo;
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public ExonExtractor(@NotNull final GeneChecker geneChecker, @NotNull final MutationTypeFilterAlgo mutationTypeFilterAlgo,
            @NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.geneChecker = geneChecker;
        this.mutationTypeFilterAlgo = mutationTypeFilterAlgo;
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @Nullable
    public List<ExonAnnotation> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (EXON_EVENTS.contains(type) && geneChecker.isValidGene(gene)) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(gene);
            assert canonicalTranscript != null;

            if (transcriptId == null || transcriptId.equals(canonicalTranscript.transcriptID())) {
                List<Integer> exonIndices = extractExonIndices(event);
                if (exonIndices == null) {
                    LOGGER.warn("Could not extract exon indices from '{}'", event);
                    return null;
                }

                MutationTypeFilter mutationTypeFilter = mutationTypeFilterAlgo.determine(gene, event);

                List<ExonAnnotation> annotations = Lists.newArrayList();
                for (int exonIndex : exonIndices) {
                    ExonAnnotation annotation = determineExonAnnotation(gene,
                            canonicalTranscript,
                            exonIndex,
                            mutationTypeFilter,
                            canonicalTranscript.transcriptID());
                    if (annotation != null) {
                        annotations.add(annotation);
                    } else {
                        LOGGER.warn("Could not determine exon annotation for exon index {} on transcript '{}' on '{}'",
                                exonIndex,
                                canonicalTranscript.transcriptID(),
                                gene);
                    }
                }
                return !annotations.isEmpty() ? annotations : null;
            } else {
                LOGGER.warn("Transcript IDs not equal for provided transcript '{}' and HMF canonical transcript '{}' for {} ",
                        transcriptId,
                        canonicalTranscript.transcriptID(),
                        event);
            }
        }

        return null;
    }

    @Nullable
    @VisibleForTesting
    static List<Integer> extractExonIndices(@NotNull String event) {
        List<Integer> exons = Lists.newArrayList();
        if (event.contains(" or ") || event.contains(" & ")) {
            exons = extractMultipleExonIndices(event);
        } else if (event.contains("-")) {
            exons = extractContinuousRangeOfExonIndices(event);
        } else {
            String[] words = event.split(" ");
            for (String word : words) {
                if (isInteger(word)) {
                    exons.add(Integer.valueOf(word));
                }
            }
        }

        return !exons.isEmpty() ? exons : null;
    }

    @NotNull
    private static List<Integer> extractMultipleExonIndices(@NotNull String event) {
        List<Integer> exonIndices = Lists.newArrayList();
        String[] words = event.replace(" or ", ",").replace(" & ", ",").replace(")", "").split(" ");
        for (String word : words) {
            if (word.contains(",")) {
                String[] exons = word.split(",");
                for (String exon : exons) {
                    exonIndices.add(Integer.valueOf(exon));
                }
            }
        }
        return exonIndices;
    }

    @NotNull
    private static List<Integer> extractContinuousRangeOfExonIndices(@NotNull String event) {
        List<Integer> exonIndices = Lists.newArrayList();
        String[] words = event.split(" ");
        for (String word : words) {
            if (word.contains("-")) {
                String[] splitEvents = word.split("-");
                int eventStart = Integer.parseInt(splitEvents[0]);
                int eventEnd = Integer.parseInt(splitEvents[1]);
                for (int i = eventStart; i <= eventEnd; i++) {
                    exonIndices.add(i);
                }
            }
        }
        return exonIndices;
    }

    @Nullable
    private static ExonAnnotation determineExonAnnotation(@NotNull String gene, @NotNull HmfTranscriptRegion transcript, int exonIndex,
            @NotNull MutationTypeFilter mutationTypeFilter, @NotNull String canonicalTranscriptID) {
        HmfExonRegion hmfExonRegion = transcript.exonByIndex(exonIndex);

        if (hmfExonRegion == null) {
            return null;
        }

        // Extend exonic range to include SPLICE variants.
        // First exon does not start with a splice region but we don't take this into account since it would not matter downstream anyways.
        long start = hmfExonRegion.start() - SPLICE_SIZE;
        long end = hmfExonRegion.end() + SPLICE_SIZE;

        return ImmutableExonAnnotation.builder()
                .gene(gene)
                .transcript(canonicalTranscriptID)
                .chromosome(hmfExonRegion.chromosome())
                .start(start)
                .end(end)
                .mutationType(mutationTypeFilter)
                .exonIndex(exonIndex)
                .build();
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
