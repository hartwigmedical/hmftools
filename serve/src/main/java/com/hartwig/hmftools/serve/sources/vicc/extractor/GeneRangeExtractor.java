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
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeType;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
                                    extractMutationTypeFilter(feature.name(), driverGenes, feature.geneSymbol()));
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
                            GeneRangeAnnotation annotation = determineCodonAnnotation(canonicalTranscript,
                                    extractMutationTypeFilter(feature.name(), driverGenes, feature.geneSymbol()),
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
                .mutationType(specificMutationType)
                .rangeType(GeneRangeType.EXON)
                .exonId(hmfExonRegion.exonID())
                .rangeNumber(exonNumber)
                .build();
    }

    @Nullable
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
    private static GeneRangeAnnotation determineCodonAnnotation(@NotNull HmfTranscriptRegion canonicalTranscript,
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
                    .mutationType(specificMutationType)
                    .rangeType(GeneRangeType.CODON)
                    .rangeNumber(codonNumber)
                    .build();
        } else if (genomeRegions != null && genomeRegions.size() == 2) {
            String chromosome = genomeRegions.get(0).chromosome().equals(genomeRegions.get(1).chromosome())
                    ? genomeRegions.get(0).chromosome()
                    : Strings.EMPTY;
            if (chromosome.isEmpty()) {
                LOGGER.warn("Chromosome positions are not equals!");
            }

            long start = canonicalTranscript.strand() == Strand.FORWARD ? genomeRegions.get(0).start() : genomeRegions.get(1).start();
            long end = canonicalTranscript.strand() == Strand.FORWARD ? genomeRegions.get(1).start() : genomeRegions.get(0).end();

            return ImmutableGeneRangeAnnotation.builder()
                    .gene(geneSymbol)
                    .chromosome(chromosome)
                    .start(start)
                    .end(end)
                    .mutationType(specificMutationType)
                    .rangeType(GeneRangeType.CODON)
                    .rangeNumber(codonNumber)
                    .build();
        } else {
            LOGGER.warn("None genomic positions of the codon {} on gene {} could be extracted", codonNumber, geneSymbol);
        }
        return null;
    }

    @VisibleForTesting
    @NotNull
    static MutationTypeFilter extractMutationTypeFilter(@NotNull String featureName, @NotNull List<DriverGene> driverGenes,
            @NotNull String gene) {
        String featureEvent = featureName.toLowerCase();
        String extractSpecificInfoOfEvent = featureEvent.substring(featureEvent.lastIndexOf(" ") + 1);
        MutationTypeFilter filter;
        if (featureEvent.contains("skipping mutation") || featureEvent.contains("splice site insertion")) {
            filter = MutationTypeFilter.SPLICE;
        } else if (extractSpecificInfoOfEvent.equals("deletions") || extractSpecificInfoOfEvent.equals("deletion") || featureEvent.contains(
                "partial deletion of exons")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_DELETION;
        } else if (extractSpecificInfoOfEvent.equals("insertions") || extractSpecificInfoOfEvent.equals("insertion")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_INSERTION;
        } else if (extractSpecificInfoOfEvent.equals("deletion/insertion") || extractSpecificInfoOfEvent.equals("insertions/deletions")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_ANY;
        } else if (extractSpecificInfoOfEvent.equals("frameshift")) {
            filter = MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        } else {
            filter = MutationTypeFilter.UNKNOWN;
        }

        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    if (filter == MutationTypeFilter.UNKNOWN) {
                        return MutationTypeFilter.MISSENSE_ANY;
                    } else {
                        return filter;
                    }
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    if (filter == MutationTypeFilter.UNKNOWN) {
                        return MutationTypeFilter.ANY;

                    } else {
                        return filter;
                    }
                }
            }
        }

        return filter;
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
