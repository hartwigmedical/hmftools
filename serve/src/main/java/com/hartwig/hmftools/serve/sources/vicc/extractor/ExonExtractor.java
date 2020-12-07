package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ExonExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ExonExtractor.class);

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public ExonExtractor(@NotNull final GeneChecker geneChecker, @NotNull final List<DriverGene> driverGenes,
            @NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.geneChecker = geneChecker;
        this.driverGenes = driverGenes;
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public Map<Feature, List<ExonAnnotation>> extract(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<ExonAnnotation>> exonsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.EXON || feature.type() == MutationType.FUSION_PAIR_AND_EXON) {
                if (geneChecker.isValidGene(feature.geneSymbol())) {
                    HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());
                    assert canonicalTranscript != null;

                    String transcriptIdVicc = viccEntry.transcriptId();
                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        List<Integer> exonNumbers = extractExonNumbers(feature.name());
                        List<ExonAnnotation> annotations = Lists.newArrayList();
                        for (int exonNumber : exonNumbers) {
                            ExonAnnotation annotation = determineExonAnnotation(feature.geneSymbol(),
                                    canonicalTranscript,
                                    exonNumber,
                                    MutationTypeFilterExtraction.extract(feature.name(),
                                            driverGenes,
                                            feature.geneSymbol()));
                            if (annotation != null) {
                                annotations.add(annotation);
                            }
                        }
                        exonsPerFeature.put(feature, annotations);
                    } else {
                        LOGGER.warn("Transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            }
        }

        return exonsPerFeature;
    }

    @NotNull
    @VisibleForTesting
    static List<Integer> extractExonNumbers(@NotNull String featureName) {
        List<Integer> exons = Lists.newArrayList();
        if (featureName.contains(" or ") || featureName.contains(" & ")) {
            exons = extractMultipleExonNumbers(featureName);
        } else if (featureName.contains("-")) {
            exons = extractListOfExonNumbers(featureName);
        } else {
            String[] words = featureName.split(" ");
            for (String word : words) {
                if (isInteger(word)) {
                    exons.add(Integer.valueOf(word));
                }
            }
        }

        if (exons.isEmpty()) {
            LOGGER.warn("Could not extract exon numbers from '{}'", featureName);
        }
        return exons;
    }

    @NotNull
    private static List<Integer> extractMultipleExonNumbers(@NotNull String featureName) {
        List<Integer> exonNumbers = Lists.newArrayList();
        String[] words = featureName.replace(" or ", ",").replace(" & ", ",").replace(")", "").split(" ");
        for (String word : words) {
            if (word.contains(",")) {
                String[] exons = word.split(",");
                for (String exon : exons) {
                    exonNumbers.add(Integer.valueOf(exon));
                }
            }
        }
        return exonNumbers;
    }

    @NotNull
    private static List<Integer> extractListOfExonNumbers(@NotNull String featureName) {
        List<Integer> exonNumbers = Lists.newArrayList();
        String[] words = featureName.split(" ");
        for (String word : words) {
            if (word.contains("-")) {
                String[] splitEvents = word.split("-");
                int eventStart = Integer.parseInt(splitEvents[0]);
                int eventEnd = Integer.parseInt(splitEvents[1]);
                for (int i = eventStart; i <= eventEnd; i++) {
                    exonNumbers.add(i);
                }
            }
        }
        return exonNumbers;
    }

    @Nullable
    private static ExonAnnotation determineExonAnnotation(@NotNull String gene, @NotNull HmfTranscriptRegion transcript, int exonIndex,
            @NotNull MutationTypeFilter specificMutationType) {
        HmfExonRegion hmfExonRegion = transcript.exonByIndex(exonIndex);

        if (hmfExonRegion == null) {
            LOGGER.warn("Could not resolve exon index {} from transcript {}", exonIndex, transcript.transcriptID());
        }

        // Extend exonic range by 5 to include SPLICE variants.
        long start = hmfExonRegion.start() - 5;
        long end = hmfExonRegion.end() + 5;

        return ImmutableExonAnnotation.builder()
                .chromosome(hmfExonRegion.chromosome())
                .start(start)
                .end(end)
                .gene(gene)
                .mutationType(specificMutationType)
                .exonEnsemblId(hmfExonRegion.exonID())
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
