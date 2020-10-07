package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class GeneRangeExtractor {

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public GeneRangeExtractor(@NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public Map<Feature, List<GeneRangeAnnotation>> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        List<GeneRangeAnnotation> geneRangeAnnotation = Lists.newArrayList();
        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());
            FeatureType featureType = feature.type();

            if (featureType == FeatureType.GENE_RANGE_EXON) {
                String transcriptIdVicc = viccEntry.transcriptId();

                if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
                    if (feature.name().contains(",")) {
                        String[] exons = feature.name()
                                .substring((feature.name().toLowerCase().indexOf("exon")))
                                .replace(" or ", ",")
                                .replace("exon ", "")
                                .split(",");
                        for (String exon : exons) {
                            int exonNumberList = Integer.parseInt(exon) - 1; // HmfExonRegion start with count 0 so exonNumber is one below
                            geneRangeAnnotation.add(extractExonGenomicPositions(feature, canonicalTranscript, exonNumberList));
                        }
                        geneRangesPerFeature.put(feature, geneRangeAnnotation);

                    } else if (feature.name().contains("or")) {
                        String[] exons = feature.name()
                                .substring((feature.name().toLowerCase().indexOf("exon")))
                                .replace(" or ", ",")
                                .replace("exon ", "")
                                .split(",");
                        for (String exon : exons) {
                            int exonNumberList = Integer.parseInt(exon) - 1; // HmfExonRegion start with count 0 so exonNumber is one below
                            geneRangeAnnotation.add(extractExonGenomicPositions(feature, canonicalTranscript, exonNumberList));

                        }
                        geneRangesPerFeature.put(feature, geneRangeAnnotation);
                    } else if (feature.description().equals("NPM1 EXON 12 MUTATION")) {
                        //Skipping because transcript has 11 exons and not 12 on the canonical transcript
                        //TODO how to solve this event?
                        //  LOGGER.warn("Skipped future for determine genomic positions of exon range '{}'", feature);
                    } else if (feature.name().contains("-")) {
                        String exons = feature.proteinAnnotation();
                        List<HmfExonRegion> exonRegions = canonicalTranscript.exome();

                        if (exons.equals("mutation")) {
                            exons = feature.name().substring((feature.name().toLowerCase().indexOf("exon"))).replace("exon ", "");
                        } else if (exons.equals("EX")) {
                            exons = feature.name().split(" ")[2].split(" ")[0];
                        }
                        int startExon =
                                Integer.parseInt(exons.split("-")[0]) - 1; // HmfExonRegion start with count 0 so exonNumber is one below
                        int endExon =
                                Integer.parseInt(exons.split("-")[1]) - 1; // HmfExonRegion start with count 0 so exonNumber is one below

                        HmfExonRegion hmfExonRegionStart = exonRegions.get(startExon);
                        HmfExonRegion hmfExonRegionEnd = exonRegions.get(endExon);

                        long start = hmfExonRegionStart.start();
                        long end = hmfExonRegionEnd.end();
                        String chromosome = hmfExonRegionStart.chromosome();

                        geneRangeAnnotation.add(ImmutableGeneRangeAnnotation.builder()
                                .gene(feature.geneSymbol())
                                .start(start)
                                .end(end)
                                .chromosome(chromosome)
                                .mutationType(MutationTypeFilter.ANY)
                                .build());
                        geneRangesPerFeature.put(feature, geneRangeAnnotation);
                    } else {

                        String exonNumber = feature.proteinAnnotation();

                        if (exonNumber.equals("mutation")) {
                            exonNumber = feature.name().substring((feature.name().toLowerCase().indexOf("exon"))).replace("exon ", "");
                        } else if (exonNumber.equals("exon")) {
                            //exon ...insertions/deletions. Determine of this is a range
                            exonNumber = feature.name()
                                    .substring((feature.name().toLowerCase().indexOf("exon")))
                                    .replace("exon ", "")
                                    .replace(" deletions", "")
                                    .replace(" insertions", "");
                        } else if (exonNumber.equals("proximal")) {
                            //check what this means
                            exonNumber = feature.name().substring((feature.name().toLowerCase().indexOf("exon"))).replace("exon ", "");
                        }
                        int exonNumberList =
                                Integer.parseInt(exonNumber) - 1; // HmfExonRegion start with count 0 so exonNumber is one below

                        geneRangeAnnotation.add(extractExonGenomicPositions(feature, canonicalTranscript, exonNumberList));

                        geneRangesPerFeature.put(feature, geneRangeAnnotation);
                    }
                } else {
                    //                    LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
                    //                            transcriptIdVicc,
                    //                            canonicalTranscript.transcriptID(),
                    //                            feature);
                }

            } else if (featureType == FeatureType.GENE_RANGE_CODON) {
                String proteinAnnotation = feature.proteinAnnotation();
                if (!proteinAnnotation.equals("T148HFSX9") && !proteinAnnotation.equals("L485_P490")) {
                    geneRangesPerFeature = determineRanges(viccEntry,
                            feature,
                            proteinAnnotation,
                            geneRangeAnnotation,
                            geneRangesPerFeature,
                            canonicalTranscript);
                }
            } else if (featureType == FeatureType.GENE_RANGE_CODON) {
                String proteinAnnotation = feature.proteinAnnotation();
                if (!proteinAnnotation.equals("T148HFSX9") && !proteinAnnotation.equals("L485_P490")) {
                    geneRangesPerFeature = determineRanges(viccEntry,
                            feature,
                            proteinAnnotation,
                            geneRangeAnnotation,
                            geneRangesPerFeature,
                            canonicalTranscript);
                }
            } else if (featureType == FeatureType.GENE_RANGE_CODON) {
                String proteinAnnotation = feature.proteinAnnotation();
                if (!proteinAnnotation.equals("T148HFSX9") && !proteinAnnotation.equals("L485_P490")) {
                    //                    geneRangesPerFeature = determineRanges(viccEntry,
                    //                            feature,
                    //                            feature.proteinAnnotation(),
                    //                            geneRangeAnnotation,
                    //                            geneRangesPerFeature,
                    //                            canonicalTranscript);

                    geneRangeAnnotation.add(ImmutableGeneRangeAnnotation.builder()
                            .gene(Strings.EMPTY)
                            .start(0)
                            .end(0)
                            .chromosome(Strings.EMPTY)
                            .mutationType(MutationTypeFilter.ANY)
                            .build());
                    geneRangesPerFeature.put(feature, geneRangeAnnotation);
                }
            } else if (featureType == FeatureType.GENE_RANGE_CODON) {
                String proteinAnnotation = feature.proteinAnnotation();
                if (!proteinAnnotation.equals("T148HFSX9") && !proteinAnnotation.equals("L485_P490")) {
                    //                    geneRangesPerFeature = determineRanges(viccEntry,
                    //                            feature,
                    //                            feature.proteinAnnotation(),
                    //                            geneRangeAnnotation,
                    //                            geneRangesPerFeature,
                    //                            canonicalTranscript);

                    geneRangeAnnotation.add(ImmutableGeneRangeAnnotation.builder()
                            .gene(Strings.EMPTY)
                            .start(0)
                            .end(0)
                            .chromosome(Strings.EMPTY)
                            .mutationType(MutationTypeFilter.ANY)
                            .build());
                    geneRangesPerFeature.put(feature, geneRangeAnnotation);
                }
            } else if (featureType == FeatureType.GENE_RANGE_CODON) {
                String proteinAnnotation = feature.proteinAnnotation();
                int start = Integer.parseInt(proteinAnnotation.split("_")[0].replaceAll("\\D+", ""));
                int end = Integer.parseInt(proteinAnnotation.split("_")[1].replaceAll("\\D+", ""));
                int count = end - start;

                if (count > 12) {
                    geneRangeAnnotation.add(ImmutableGeneRangeAnnotation.builder()
                            .gene(Strings.EMPTY)
                            .start(0)
                            .end(0)
                            .chromosome(Strings.EMPTY)
                            .mutationType(MutationTypeFilter.ANY)
                            .build());
                    geneRangesPerFeature.put(feature, geneRangeAnnotation);
                }
            }
        }
        return geneRangesPerFeature;
    }

    @NotNull
    private static GeneRangeAnnotation extractExonGenomicPositions(@NotNull Feature feature,
            @NotNull HmfTranscriptRegion canonicalTranscript, int exonNumberList) {
        List<HmfExonRegion> exonRegions = canonicalTranscript.exome();
        HmfExonRegion hmfExonRegion = exonRegions.get(exonNumberList);
        long start = hmfExonRegion.start();
        long end = hmfExonRegion.end();
        String chromosome = hmfExonRegion.chromosome();

        return ImmutableGeneRangeAnnotation.builder()
                .gene(feature.geneSymbol())
                .start(start)
                .end(end)
                .chromosome(chromosome)
                .mutationType(MutationTypeFilter.ANY)
                .build();
    }

    private static Map<Feature, List<GeneRangeAnnotation>> determineRanges(@NotNull ViccEntry viccEntry, @NotNull Feature feature,
            @NotNull String proteinAnnotation, @NotNull List<GeneRangeAnnotation> geneRangeAnnotation,
            @NotNull Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature, @NotNull HmfTranscriptRegion canonicalTranscript) {
        String transcriptIdVicc = viccEntry.transcriptId();

        //        if (transcriptIdVicc == null || transcriptIdVicc.equals(canonicalTranscript.transcriptID())) {
        //            String geneSymbol = feature.geneSymbol();
        //            int codonNumber = Integer.valueOf(proteinAnnotation.replaceAll("\\D+", ""));
        //            List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonNumber);
        //            if (genomeRegions.size() == 1) {
        //                long start = genomeRegions.get(0).start();
        //                long end = genomeRegions.get(0).end();
        //                String chromosome = genomeRegions.get(0).chromosome();
        //
        //                geneRangeAnnotation.add(ImmutableGeneRangeAnnotation.builder()
        //                        .gene(geneSymbol)
        //                        .start(start)
        //                        .end(end)
        //                        .chromosome(chromosome)
        //                        .event(feature.name())
        //                        .build());
        //                geneRangesPerFeature.put(feature, geneRangeAnnotation);
        //
        //            } else {
        //                LOGGER.warn("Multiple genomic regions known for event {}", feature);
        //            }
        //        } else {
        //            //                    LOGGER.warn("transcript IDs not equal for transcript VICC {} and HMF {} for {} ",
        //            //                            transcriptIdVicc,
        //            //                            canonicalTranscript.transcriptID(),
        //            //                            feature);
        //        }
        return geneRangesPerFeature;
    }

    private static boolean isValidSingleCodonRange(@NotNull String feature) {

        // Features are expected to look something like V600 (1 char - N digits)
        if (feature.length() < 3) {
            return false;
        }

        if (!Character.isLetter(feature.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(feature.charAt(1))) {
            return false;
        }

        if (feature.contains("*")) {
            return false;
        }

        if (feature.contains("/")) {
            return false;
        }

        if (feature.contains("fs")) {
            return false;
        }

        return Character.isDigit(feature.substring(feature.length() - 1).charAt(0));
    }
}
