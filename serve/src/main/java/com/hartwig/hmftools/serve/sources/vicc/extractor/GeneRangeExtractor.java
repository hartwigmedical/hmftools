package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneRangeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneRangeExtractor.class);

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final GeneChecker geneChecker;

    public GeneRangeExtractor(@NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap,
            @NotNull final List<DriverGene> driverGenes, @NotNull final GeneChecker geneChecker) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
        this.driverGenes = driverGenes;
        this.geneChecker = geneChecker;
    }

    @NotNull
    public Map<Feature, List<GeneRangeAnnotation>> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

            if (feature.type() == MutationType.EXON || feature.type() == MutationType.FUSION_PAIR_AND_EXON) {
                if (geneChecker.isValidGene(feature.geneSymbol(), canonicalTranscript, feature.name(), null)) {
                    String transcriptIdVicc = viccEntry.transcriptId();
                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        List<Integer> exonNumbers = extractExonNumbers(feature.name());
                        List<GeneRangeAnnotation> annotations = Lists.newArrayList();
                        for (int exonNumber : exonNumbers) {
                            GeneRangeAnnotation annotation = determineExonAnnotation(feature.geneSymbol(),
                                    canonicalTranscript,
                                    exonNumber,
                                    extractSpecificMutationTypeFilter(feature.name()));
                            if (annotation != null) {
                                annotations.add(annotation);
                            }
                        }
                        geneRangesPerFeature.put(feature, annotations);
                    } else {
                        LOGGER.warn("Transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            } else if (feature.type() == MutationType.CODON) {
                if (geneChecker.isValidGene(feature.geneSymbol(), canonicalTranscript, feature.name(), null)) {
                    String transcriptIdVicc = viccEntry.transcriptId();

                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        Integer codonNumber = extractCodonNumber(feature.name());
                        if (codonNumber != null) {
                            List<GeneRangeAnnotation> annotations = Lists.newArrayList();
                            String geneSymbol = feature.geneSymbol();
                            GeneRangeAnnotation annotation = determineCodonAnnotation(feature,
                                    canonicalTranscript,
                                    driverGenes,
                                    extractSpecificMutationTypeFilter(feature.name()),
                                    codonNumber,
                                    geneSymbol);
                            if (annotation != null) {
                                annotations.add(annotation);
                            }
                            geneRangesPerFeature.put(feature, annotations);
                        }

                    } else {
                        LOGGER.warn("Transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            }
        }

        return geneRangesPerFeature;
    }

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
            LOGGER.warn("No exon number is extracted from event '{}'", featureName);
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

    @VisibleForTesting
    static Integer extractCodonNumber(@NotNull String featureName) {
        String codonPart;
        if (featureName.contains(" ")) {
            codonPart = featureName.split(" ")[1];
        } else {
            codonPart = featureName;
        }
        codonPart = codonPart.replaceAll("\\D+", "");
        if (isInteger(codonPart)) {
            return Integer.parseInt(codonPart);
        } else {
            LOGGER.warn("Could not extract codon number from '{}'", featureName);
            return null;
        }
    }

    @Nullable
    private static GeneRangeAnnotation determineExonAnnotation(@NotNull String gene, @NotNull HmfTranscriptRegion transcript,
            int exonNumber, @NotNull MutationTypeFilter specificMutationType) {
        HmfExonRegion hmfExonRegion = transcript.exonByIndex(exonNumber);

        if (hmfExonRegion == null) {
            LOGGER.warn("Could not resolve exon {} from transcript {}", exonNumber, transcript.transcriptID());
        }

        // Extend exonic range by 5 to include SPLICE variants.
        long start = hmfExonRegion.start() - 5;
        long end = hmfExonRegion.end() + 5;

        return ImmutableGeneRangeAnnotation.builder()
                .gene(gene)
                .chromosome(hmfExonRegion.chromosome())
                .start(start)
                .end(end)
                .rangeInfo(exonNumber)
                .exonId(hmfExonRegion.exonID())
                .mutationType(specificMutationType)
                .build();
    }

    @Nullable
    private static ImmutableGeneRangeAnnotation determineCodonAnnotation(@NotNull Feature feature,
            @NotNull HmfTranscriptRegion canonicalTranscript, @NotNull List<DriverGene> driverGenes,
            @NotNull MutationTypeFilter specificMutationType, int codonNumber, @NotNull String geneSymbol) {
        List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonNumber);
        if (genomeRegions != null && genomeRegions.size() == 1) {
            String chromosome = genomeRegions.get(0).chromosome();
            long start = genomeRegions.get(0).start();
            long end = genomeRegions.get(0).end();

            return ImmutableGeneRangeAnnotation.builder()
                    .gene(geneSymbol)
                    .chromosome(chromosome)
                    .start(start)
                    .end(end)
                    .rangeInfo(codonNumber)
                    .mutationType(extractMutationFilter(driverGenes, geneSymbol, specificMutationType, feature))
                    .build();
        }

        return null;
    }

    @VisibleForTesting
    @NotNull
    static MutationTypeFilter extractMutationFilter(@NotNull List<DriverGene> driverGenes, @NotNull String gene,
            @NotNull MutationTypeFilter specificMutationType, @NotNull Feature feature) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    if (specificMutationType == MutationTypeFilter.UNKNOWN) {
                        return MutationTypeFilter.MISSENSE_ANY;
                    } else {
                        return specificMutationType;
                    }
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    if (specificMutationType == MutationTypeFilter.UNKNOWN) {
                        //TODO which inactivation event for TSG?
                        if (feature.biomarkerType() != null && feature.biomarkerType().contains("")) {
                            return MutationTypeFilter.MISSENSE_ANY;
                        } else {
                            return MutationTypeFilter.ANY;
                        }
                    } else {
                        return specificMutationType;
                    }
                }
            }
        }
        return MutationTypeFilter.UNKNOWN;
    }

    @VisibleForTesting
    @NotNull
    static MutationTypeFilter extractSpecificMutationTypeFilter(@NotNull String featureName) {
        String featureEvent = featureName.toLowerCase();
        String extractSpecificInfoOfEvent = featureEvent.substring(featureEvent.lastIndexOf(" ") + 1);
        if (featureEvent.contains("skipping mutation") || featureEvent.contains("splice site insertion")) {
            return MutationTypeFilter.SPLICE;
        } else if (extractSpecificInfoOfEvent.equals("deletions") || extractSpecificInfoOfEvent.equals("deletion") || featureEvent.contains(
                "partial deletion of exons")) {
            return MutationTypeFilter.MISSENSE_INFRAME_DELETION;
        } else if (extractSpecificInfoOfEvent.equals("insertions") || extractSpecificInfoOfEvent.equals("insertion")) {
            return MutationTypeFilter.MISSENSE_INFRAME_INSERTION;
        } else if (extractSpecificInfoOfEvent.equals("deletion/insertion") || extractSpecificInfoOfEvent.equals("insertions/deletions")) {
            return MutationTypeFilter.MISSENSE_INFRAME_ANY;
        } else if (extractSpecificInfoOfEvent.equals("frameshift")) {
            return MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        }

        return MutationTypeFilter.UNKNOWN;
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
