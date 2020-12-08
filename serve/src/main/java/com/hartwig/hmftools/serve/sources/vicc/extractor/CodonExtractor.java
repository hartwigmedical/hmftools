package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.GeneChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CodonExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CodonExtractor.class);

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public CodonExtractor(@NotNull final GeneChecker geneChecker, @NotNull final List<DriverGene> driverGenes,
            @NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.geneChecker = geneChecker;
        this.driverGenes = driverGenes;
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public List<CodonAnnotation> extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        if (type == EventType.CODON && geneChecker.isValidGene(gene)) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(gene);
            assert canonicalTranscript != null;

            if (transcriptId == null || transcriptId.equals(canonicalTranscript.transcriptID())) {
                Integer codonNumber = extractCodonNumber(event);
                if (codonNumber != null) {
                    List<CodonAnnotation> annotations = Lists.newArrayList();
                    CodonAnnotation annotation = determineCodonAnnotation(canonicalTranscript,
                            MutationTypeFilterExtraction.extract(event, driverGenes, gene),
                            codonNumber,
                            gene);
                    if (annotation != null) {
                        annotations.add(annotation);
                    }
                    return annotations;
                }
            } else {
                LOGGER.warn("Transcript IDs not equal for provided transcript '{}' and HMF canonical transcript '{}' for {} ",
                        transcriptId,
                        canonicalTranscript.transcriptID(),
                        event);
            }
        }

        return Lists.newArrayList();
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
    private static CodonAnnotation determineCodonAnnotation(@NotNull HmfTranscriptRegion canonicalTranscript,
            @NotNull MutationTypeFilter specificMutationType, int codonIndex, @NotNull String geneSymbol) {
        List<GenomeRegion> genomeRegions = canonicalTranscript.codonByIndex(codonIndex);

        if (genomeRegions != null && genomeRegions.size() == 1) {
            String chromosome = genomeRegions.get(0).chromosome();
            long start = genomeRegions.get(0).start();
            long end = genomeRegions.get(0).end();

            return ImmutableCodonAnnotation.builder()
                    .chromosome(chromosome)
                    .start(start)
                    .end(end)
                    .gene(geneSymbol)
                    .mutationType(specificMutationType)
                    .codonIndex(codonIndex)
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

            return ImmutableCodonAnnotation.builder()
                    .chromosome(chromosome)
                    .start(start)
                    .end(end)
                    .gene(geneSymbol)
                    .mutationType(specificMutationType)
                    .codonIndex(codonIndex)
                    .build();
        } else {
            LOGGER.warn("None genomic positions of the codon {} on gene {} could be extracted", codonIndex, geneSymbol);
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
