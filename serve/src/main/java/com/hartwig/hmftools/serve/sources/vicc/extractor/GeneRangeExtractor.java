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
import com.hartwig.hmftools.serve.sources.vicc.check.CheckGenes;
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
        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

            if (feature.type() == MutationType.EXON || feature.type() == MutationType.FUSION_PAIR_AND_EXON) {
                if (canonicalTranscript == null) {
                    CheckGenes.checkGensInPanel(feature.geneSymbol(), feature.name());
                } else {
                    String transcriptIdVicc = viccEntry.transcriptId();
                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        List<Integer> exonNumbers = extractExonNumbers(feature.name());
                        List<GeneRangeAnnotation> annotations = Lists.newArrayList();
                        for (int exonNumber : exonNumbers) {
                            annotations.add(extractGeneRangeAnnotationForFeature(exonNumber,
                                    feature,
                                    canonicalTranscript,
                                    driverGenes,
                                    extractSpecificMutationTypeFilter(feature)));
                        }
                        geneRangesPerFeature.put(feature, annotations);
                    } else {
                        LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                                transcriptIdVicc,
                                canonicalTranscript.transcriptID(),
                                feature);
                    }
                }
            } else if (feature.type() == MutationType.CODON) {
                if (canonicalTranscript == null) {
                    CheckGenes.checkGensInPanel(feature.geneSymbol(), feature.name());
                } else {
                    String transcriptIdVicc = viccEntry.transcriptId();

                    if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                        Integer codonNumber = extractCodonNumber(feature.name());
                        List<GeneRangeAnnotation> annotations = Lists.newArrayList();
                        if (codonNumber != null) {
                            String geneSymbol = feature.geneSymbol();
                            annotations.add(determineRanges(feature,
                                    canonicalTranscript,
                                    driverGenes,
                                    extractSpecificMutationTypeFilter(feature),
                                    codonNumber,
                                    geneSymbol));
                            geneRangesPerFeature.put(feature, annotations);
                        }

                    } else {
                        LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
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
    private static ImmutableGeneRangeAnnotation determineRanges(@NotNull Feature feature, @NotNull HmfTranscriptRegion canonicalTranscript,
            @NotNull List<DriverGene> driverGenes, @NotNull MutationTypeFilter specificMutationType, int codonNumber,
            @NotNull String geneSymbol) {
        ImmutableGeneRangeAnnotation.Builder geneRangeAnnotation = ImmutableGeneRangeAnnotation.builder();

        List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonNumber);
        if (genomeRegions != null && genomeRegions.size() == 1) {
            long start = genomeRegions.get(0).start();
            long end = genomeRegions.get(0).end();
            String chromosome = genomeRegions.get(0).chromosome();

            geneRangeAnnotation
                    .gene(geneSymbol)
                    .start(start)
                    .end(end)
                    .chromosome(chromosome)
                    .rangeInfo(codonNumber)
                    .mutationType(extractMutationFilter(driverGenes, geneSymbol, specificMutationType, feature));
        }
        return geneRangeAnnotation.build();
    }

    @VisibleForTesting
    @NotNull
    public static MutationTypeFilter extractMutationFilter(@NotNull List<DriverGene> driverGenes, @NotNull String gene,
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
        CheckGenes.checkGensInPanel(gene, feature.name());
        return MutationTypeFilter.UNKNOWN;
    }

    @VisibleForTesting
    @NotNull
    public static MutationTypeFilter extractSpecificMutationTypeFilter(@NotNull Feature feature) {
        String featureEvent = feature.name().toLowerCase();
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

    @NotNull
    private static GeneRangeAnnotation extractGeneRangeAnnotationForFeature(int exonNumber, @NotNull Feature feature,
            @NotNull HmfTranscriptRegion canonicalTranscript, @NotNull List<DriverGene> driverGenes,
            @NotNull MutationTypeFilter specificMutationType) {
        int exonNumberList = exonNumber - 1; // HmfExonRegion start with count 0 so exonNumber is one below
        return extractExonGenomicPositions(feature, canonicalTranscript, exonNumberList, driverGenes, exonNumber, specificMutationType);
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
