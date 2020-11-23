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
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneRangeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneRangeExtractor.class);

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;
    @NotNull
    private final List<DriverGene> driverGenes;

    public GeneRangeExtractor(@NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap,
            @NotNull final List<DriverGene> driverGenes) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
        this.driverGenes = driverGenes;
    }

    @NotNull
    public Map<Feature, List<GeneRangeAnnotation>> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        List<GeneRangeAnnotation> geneRangeAnnotation = Lists.newArrayList();
        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

            MutationType mutationType = feature.type();

            if (mutationType == MutationType.FUSION_PAIR_AND_GENE_RANGE_EXON) {
                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene {} in HMF gene panel. Skipping fusion pair and gene range extraction for range!",
                            feature.geneSymbol());
                } else {
                    String transcriptIdVicc = viccEntry.transcriptId();
                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {

                        List<Integer> exonNumbers = extractExonNumbers(feature.name());
                        for (int exonNumber : exonNumbers) {
                            extractGeneRangesPerFeature(exonNumber,
                                    feature,
                                    canonicalTranscript,
                                    driverGenes,
                                    geneRangeAnnotation,
                                    geneRangesPerFeature,
                                    extractSpecificMutationTypeFilter(feature));
                        }
                    } else {
                        LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            }
            if (mutationType == MutationType.GENE_RANGE_EXON) {
                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene {} in HMF gene panel. Skipping gene range exon extraction!", feature.geneSymbol());
                } else {
                    String transcriptIdVicc = viccEntry.transcriptId();
                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        List<Integer> exonNumbers = extractExonNumbers(feature.name());
                        for (int exonNumber : exonNumbers) {
                            extractGeneRangesPerFeature(exonNumber,
                                    feature,
                                    canonicalTranscript,
                                    driverGenes,
                                    geneRangeAnnotation,
                                    geneRangesPerFeature,
                                    extractSpecificMutationTypeFilter(feature));
                        }
                    } else {
                        LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            } else if (mutationType == MutationType.GENE_RANGE_CODON) {
                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene {} in HMF gene panel. Skipping gene range codon extraction!", feature.geneSymbol());
                } else {
                    Integer codonNumber = extractCodonNumber(feature.name());
                    if (codonNumber != null) {
                        String geneSymbol = feature.geneSymbol();
                        geneRangesPerFeature = determineRanges(viccEntry,
                                feature,
                                geneRangeAnnotation,
                                geneRangesPerFeature,
                                canonicalTranscript,
                                driverGenes,
                                extractSpecificMutationTypeFilter(feature),
                                codonNumber,
                                geneSymbol);
                        geneRangesPerFeature.put(feature, geneRangeAnnotation);
                    }
                }
            }
        }

        return geneRangesPerFeature;
    }

    @VisibleForTesting
    static List<Integer> extractExonNumbers(@NotNull String featureName) {
        List<Integer> codons = Lists.newArrayList();
        if (featureName.contains("or") && featureName.contains(",")) {
            String formattedFeatureName = featureName.replace(" or ", ",");
            String[] splitEvents = formattedFeatureName.split(" ")[4].split(",");
            for (String events : splitEvents) {
                codons.add(Integer.valueOf(events));
            }
        } else if (featureName.contains("or")) {
            String formattedFeatureName = featureName.replace(" or ", ",");
            String[] splitEvents = formattedFeatureName.split(" ")[4].split(",");
            for (String events : splitEvents) {
                codons.add(Integer.valueOf(events));
            }
        } else if (featureName.contains("-")) {
            String[] splitEvents = featureName.split("-");
            String[] eventStartExtract = splitEvents[0].split(" ");
            int eventStart = Integer.parseInt(eventStartExtract[eventStartExtract.length - 1]);
            String[] eventEndExtract = splitEvents[1].split(" ");
            int eventEnd = Integer.parseInt(eventEndExtract[0]);
            for (int i = eventStart; i <= eventEnd; i++) {
                codons.add(i);
            }
        } else if (featureName.contains("&")) {
            String formattedFeatureName = featureName.replace(" & ", ",").replace(")", "");
            String[] splitEvents = formattedFeatureName.split(" ")[5].split(",");
            for (String events : splitEvents) {
                codons.add(Integer.valueOf(events));
            }
        } else {
            String[] extraction = featureName.split(" ");
            for (String event : extraction) {
                if (event.matches("[0-9]+")) {
                    codons.add(Integer.valueOf(event));
                }
            }
        }

        if (codons.size() == 0) {
            LOGGER.warn("No exon number is extracted from event {}", featureName);
        }
        return codons;
    }

    @VisibleForTesting
    public static Integer extractCodonNumber(@NotNull String featureName) {
        if (featureName.split(" ").length == 1) {
            if (!isInteger(featureName.replaceAll("\\D+", ""))) {
                LOGGER.warn("Could not convert gene range codon {} to codon number", featureName);
                return null;
            } else {
                return Integer.parseInt(featureName.replaceAll("\\D+", ""));
            }
        } else if (featureName.split(" ").length == 2) {
            String featureNamePart = featureName.split(" ")[1];
            if (!isInteger(featureNamePart.replaceAll("\\D+", ""))) {
                LOGGER.warn("Could not convert gene range codon {} to codon number", featureName);
                return null;
            } else {
                return Integer.parseInt(featureNamePart.replaceAll("\\D+", ""));
            }
        }
        LOGGER.warn("Could not convert gene range codon {} to codon number", featureName);
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

    @NotNull
    private static Map<Feature, List<GeneRangeAnnotation>> determineRanges(@NotNull ViccEntry viccEntry, @NotNull Feature feature,
            @NotNull List<GeneRangeAnnotation> geneRangeAnnotations, @NotNull Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature,
            @NotNull HmfTranscriptRegion canonicalTranscript, @NotNull List<DriverGene> driverGenes,
            @NotNull MutationTypeFilter specificMutationType, int codonNumber, @NotNull String geneSymbol) {
        String transcriptIdVicc = viccEntry.transcriptId();

        if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {

            List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonNumber);
            if (genomeRegions != null && genomeRegions.size() == 1) {
                long start = genomeRegions.get(0).start();
                long end = genomeRegions.get(0).end();
                String chromosome = genomeRegions.get(0).chromosome();

                geneRangeAnnotations.add(ImmutableGeneRangeAnnotation.builder()
                        .gene(geneSymbol)
                        .start(start)
                        .end(end)
                        .chromosome(chromosome)
                        .rangeInfo(codonNumber)
                        .mutationType(extractMutationFilter(driverGenes, geneSymbol, specificMutationType, feature))
                        .build());
                geneRangesPerFeature.put(feature, geneRangeAnnotations);
            }
        } else {
            LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                    transcriptIdVicc,
                    canonicalTranscript.transcriptID(),
                    feature);
        }
        return geneRangesPerFeature;
    }

    @NotNull
    private static MutationTypeFilter extractMutationFilter(@NotNull List<DriverGene> driverGenes, @NotNull String gene,
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
                        if (feature.biomarkerType() != null && feature.biomarkerType().contains("inactivation")) {
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
        LOGGER.warn("Gene {} is not present in driver catalog in gene range extractor", gene);
        return MutationTypeFilter.UNKNOWN;
    }

    @NotNull
    private static MutationTypeFilter extractSpecificMutationTypeFilter(@NotNull Feature feature) {
        String featureEvent = feature.name().toLowerCase();
        String extractSpecificInfoOfEvent = featureEvent.substring(featureEvent.lastIndexOf(" ") + 1);
        if (featureEvent.contains("skipping mutation") || featureEvent.contains("splice site insertion")) {
            return MutationTypeFilter.SPLICE;
        } else if (extractSpecificInfoOfEvent.equals("deletions") || extractSpecificInfoOfEvent.equals("deletion")
                || extractSpecificInfoOfEvent.contains("Partial deletion of Exons")) {
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

    private static void extractGeneRangesPerFeature(int exonNumber, @NotNull Feature feature,
            @NotNull HmfTranscriptRegion canonicalTranscript, @NotNull List<DriverGene> driverGenes,
            @NotNull List<GeneRangeAnnotation> geneRangeAnnotation, @NotNull Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature,
            @NotNull MutationTypeFilter specificMutationType) {
        int exonNumberList = exonNumber - 1; // HmfExonRegion start with count 0 so exonNumber is one below

        geneRangeAnnotation.add(extractExonGenomicPositions(feature,
                canonicalTranscript,
                exonNumberList,
                driverGenes,
                exonNumber,
                specificMutationType));
        geneRangesPerFeature.put(feature, geneRangeAnnotation);
    }

    @NotNull
    private static GeneRangeAnnotation extractExonGenomicPositions(@NotNull Feature feature,
            @NotNull HmfTranscriptRegion canonicalTranscript, int exonNumberList, @NotNull List<DriverGene> driverGenes, int exonNumber,
            @NotNull MutationTypeFilter specificMutationType) {
        List<HmfExonRegion> exonRegions = canonicalTranscript.exome();
        HmfExonRegion hmfExonRegion = exonRegions.get(exonNumberList);
        long start = hmfExonRegion.start() - 5;
        long end = hmfExonRegion.end() + 5;
        String chromosome = hmfExonRegion.chromosome();

        return ImmutableGeneRangeAnnotation.builder()
                .gene(feature.geneSymbol())
                .start(start)
                .end(end)
                .chromosome(chromosome)
                .rangeInfo(exonNumber)
                .exonId(hmfExonRegion.exonID())
                .mutationType(extractMutationFilter(driverGenes, feature.geneSymbol(), specificMutationType, feature))
                .build();
    }
}
