package com.hartwig.hmftools.common.genome.genepanel;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.genome.region.ModifiableHmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HmfGenomeFileLoader {

    private static final Logger LOGGER = LogManager.getLogger(HmfGenomeFileLoader.class);

    private static final String FIELD_SEPARATOR = "\t";

    private static final int CHROMOSOME_COLUMN = 0;
    private static final int GENE_START_COLUMN = 1;
    private static final int GENE_END_COLUMN = 2;
    private static final int GENE_ID_COLUMN = 3;
    private static final int GENE_COLUMN = 4;
    private static final int ENTREZ_ID_COLUMN = 5;
    private static final int CHROMOSOME_BAND_COLUMN = 6;
    private static final int TRANSCRIPT_ID_COLUMN = 7;
    private static final int TRANSCRIPT_VERSION_COLUMN = 8;
    private static final int TRANSCRIPT_START_COLUMN = 9;
    private static final int TRANSCRIPT_END_COLUMN = 10;
    private static final int EXON_ID_COLUMN = 11;
    private static final int EXON_START_COLUMN = 12;
    private static final int EXON_END_COLUMN = 13;
    private static final int STRAND_COLUMN = 14;
    private static final int CODING_START_COLUMN = 15;
    private static final int CODING_END_COLUMN = 16;

    private HmfGenomeFileLoader() {
    }

    @NotNull
    public static List<HmfTranscriptRegion> fromInputStream(@NotNull InputStream genomeInputStream) {
        return fromLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static List<HmfTranscriptRegion> fromLines(@NotNull List<String> lines) {
        Map<String, HmfTranscriptRegion> geneMap = Maps.newLinkedHashMap();
        for (final String line : lines) {
            String[] values = line.split(FIELD_SEPARATOR);
            String chromosome = values[CHROMOSOME_COLUMN].trim();
            if (!HumanChromosome.contains(chromosome)) {
                LOGGER.warn("Skipping line due to unknown chromosome: {}", line);
                continue;
            }

            String gene = values[GENE_COLUMN];
            long transcriptStart = Long.parseLong(values[TRANSCRIPT_START_COLUMN].trim());
            long transcriptEnd = Long.parseLong(values[TRANSCRIPT_END_COLUMN].trim());

            if (transcriptEnd < transcriptStart) {
                LOGGER.warn("Invalid transcript region found on chromosome {}: start={}, end={}",
                        chromosome,
                        transcriptStart,
                        transcriptEnd);

            } else {
                HmfTranscriptRegion geneRegion = geneMap.computeIfAbsent(gene,
                        geneName -> createRegion(chromosome, transcriptStart, transcriptEnd, geneName, values));

                HmfExonRegion exonRegion = ImmutableHmfExonRegion.builder()
                        .chromosome(chromosome)
                        .exonRank(geneRegion.exome().size() + 1)
                        .start(Long.parseLong(values[EXON_START_COLUMN]))
                        .end(Long.parseLong(values[EXON_END_COLUMN]))
                        .build();

                geneRegion.exome().add(exonRegion);
            }
        }

        return Lists.newArrayList(geneMap.values());
    }

    @NotNull
    private static HmfTranscriptRegion createRegion(@NotNull String chromosome, long transcriptStart, long transcriptEnd,
            @NotNull String gene, @NotNull String[] values) {
        final String entrezIdString = values[ENTREZ_ID_COLUMN];

        final List<Integer> entrezIds = entrezIdString.isEmpty()
                ? Lists.newArrayList()
                : Arrays.stream(entrezIdString.split(",")).map(Integer::parseInt).collect(Collectors.toList());

        long codingStart = 0;
        long codingEnd = 0;
        if (values.length > CODING_END_COLUMN) {
            codingStart = Long.parseLong(values[CODING_START_COLUMN]);
            codingEnd = Long.parseLong(values[CODING_END_COLUMN]);
        }
        // TODO: Remove dependency on modifiable transcript region.
        return ModifiableHmfTranscriptRegion.create()
                .setChromosome(chromosome)
                .setStart(transcriptStart)
                .setEnd(transcriptEnd)
                .setTranscriptID(values[TRANSCRIPT_ID_COLUMN])
                .setTranscriptVersion(Integer.parseInt(values[TRANSCRIPT_VERSION_COLUMN]))
                .setChromosomeBand(values[CHROMOSOME_BAND_COLUMN])
                .setEntrezId(entrezIds)
                .setGene(gene)
                .setGeneID(values[GENE_ID_COLUMN])
                .setGeneStart(Long.parseLong(values[GENE_START_COLUMN]))
                .setGeneEnd(Long.parseLong(values[GENE_END_COLUMN]))
                .setStrand(Strand.valueOf(Integer.parseInt(values[STRAND_COLUMN])))
                .setCodingStart(codingStart)
                .setCodingEnd(codingEnd);
    }
}
