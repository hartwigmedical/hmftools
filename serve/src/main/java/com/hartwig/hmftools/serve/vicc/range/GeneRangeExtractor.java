package com.hartwig.hmftools.serve.vicc.range;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class GeneRangeExtractor {
    private static final Logger LOGGER = LogManager.getLogger(GeneRangeExtractor.class);

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    private static final Set<String> GENE_EXON = Sets.newHashSet("exon");
    private static final Set<String> GENE_MULTIPLE_CODONS = Sets.newHashSet("nonsense", "V600E/K");

    public GeneRangeExtractor(@NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public Map<Feature, GeneRangeAnnotation> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, GeneRangeAnnotation> geneRangesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            String event = Strings.EMPTY;
            if (feature.name().toLowerCase().contains("exon")) {
                event = "exon";
            }

            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

            if (GENE_EXON.contains(feature.name().toLowerCase()) || GENE_EXON.contains(event)) {
                //TODO: fix is ignore events
                if (!feature.name().contains(",") && !feature.name().contains("/") && !feature.name().equals("3' EXON DELETION")) {
                    String exonNumber = feature.name().substring((feature.name().toLowerCase().indexOf("exon"))).replaceAll("\\D+", "");
                    int exonNumberList = Integer.valueOf(exonNumber)-1; // HmfExonRegion start with count 0 so exonNumber is one below
                    List<HmfExonRegion> exonRegions = canonicalTranscript.exome();
                    HmfExonRegion hmfExonRegion = exonRegions.get(exonNumberList);

                    long start = hmfExonRegion.start();
                    long end = hmfExonRegion.end();
                    String chromosome = hmfExonRegion.chromosome();
                    String exonId = hmfExonRegion.exonID();

                    GeneRangeAnnotation geneRangeAnnotation = ImmutableGeneRangeAnnotation.builder()
                            .gene(feature.geneSymbol())
                            .start(start)
                            .end(end)
                            .chromosome(chromosome)
                            .event(feature.name())
                            .build();

                    geneRangesPerFeature.put(feature, geneRangeAnnotation);
                }

            } else if (GENE_MULTIPLE_CODONS.contains(feature.biomarkerType()) && feature.proteinAnnotation()
                    .substring(feature.proteinAnnotation().length() - 1)
                    .equals("X") || GENE_MULTIPLE_CODONS.contains(feature.proteinAnnotation())) {
                String geneSymbol = feature.geneSymbol();
                String proteinAnnotation = feature.proteinAnnotation();
                int codonNumber = Integer.valueOf(proteinAnnotation.replaceAll("\\D+", ""));

                List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonNumber);
                if (genomeRegions.size() == 1) {
                    long start = genomeRegions.get(0).start();
                    long end = genomeRegions.get(0).end();
                    String chromosome = genomeRegions.get(0).chromosome();

                    GeneRangeAnnotation geneRangeAnnotation = ImmutableGeneRangeAnnotation.builder()
                            .gene(geneSymbol)
                            .start(start)
                            .end(end)
                            .chromosome(chromosome)
                            .event(feature.name())
                            .build();
                    geneRangesPerFeature.put(feature, geneRangeAnnotation);

                } else {
                    LOGGER.warn("Multiple genomic regions known for event {}", feature);
                }
            }
        }
        return geneRangesPerFeature;
    }
}
