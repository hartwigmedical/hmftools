package com.hartwig.hmftools.common.region.hmfslicer;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

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
    private static final int CODING_START_COLUMN = 14;
    private static final int CODING_END_COLUMN = 15;

    private HmfGenomeFileLoader() {
    }

    @NotNull
    public static SortedSetMultimap<String, HmfGenomeRegion> fromInputStream(@NotNull final InputStream genomeInputStream) {
        return fromLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static SortedSetMultimap<String, HmfGenomeRegion> fromLines(@NotNull final List<String> lines) {
        final Map<String, ModifiableHmfGenomeRegion> geneMap = Maps.newLinkedHashMap();
        for (final String line : lines) {
            final String[] values = line.split(FIELD_SEPARATOR);
            final String chromosome = values[CHROMOSOME_COLUMN].trim();
            if (!HumanChromosome.contains(chromosome)) {
                LOGGER.warn("skipping line due to unknown chromosome: {}", line);
                continue;
            }

            final String gene = values[GENE_COLUMN];
            final long transcriptStart = Long.valueOf(values[TRANSCRIPT_START_COLUMN].trim());
            final long transcriptEnd = Long.valueOf(values[TRANSCRIPT_END_COLUMN].trim());

            if (transcriptEnd < transcriptStart) {
                LOGGER.warn("Invalid transcript region found on chromosome " + chromosome + ": start=" + transcriptStart + ", end="
                        + transcriptEnd);
            } else {
                final ModifiableHmfGenomeRegion geneRegion = geneMap.computeIfAbsent(gene,
                        geneName -> createRegion(chromosome, transcriptStart, transcriptEnd, geneName, values));

                final HmfExonRegion exonRegion = ImmutableHmfExonRegion.builder()
                        .chromosome(chromosome)
                        .exonID(values[EXON_ID_COLUMN])
                        .start(Long.valueOf(values[EXON_START_COLUMN]))
                        .end(Long.valueOf(values[EXON_END_COLUMN]))
                        .build();

                geneRegion.addExome(exonRegion);
            }
        }

        final SortedSetMultimap<String, HmfGenomeRegion> regionMap = TreeMultimap.create();
        for (ModifiableHmfGenomeRegion geneRegion : geneMap.values()) {
            regionMap.put(geneRegion.chromosome(), geneRegion);
        }

        return regionMap;
    }

    @NotNull
    private static ModifiableHmfGenomeRegion createRegion(final String chromosome, final long transcriptStart, final long transcriptEnd,
            @NotNull final String gene, @NotNull final String[] values) {
        final List<Integer> entrezIds =
                Arrays.stream(values[ENTREZ_ID_COLUMN].split(",")).map(Integer::parseInt).collect(Collectors.toList());

        long codingStart = 0;
        long codingEnd = 0;
        if (values.length > CODING_END_COLUMN) {
            codingStart = Long.valueOf(values[CODING_START_COLUMN]);
            codingEnd = Long.valueOf(values[CODING_END_COLUMN]);
        }


        return ModifiableHmfGenomeRegion.create()
                .setChromosome(chromosome)
                .setStart(transcriptStart)
                .setEnd(transcriptEnd)
                .setTranscriptID(values[TRANSCRIPT_ID_COLUMN])
                .setTranscriptVersion(Integer.valueOf(values[TRANSCRIPT_VERSION_COLUMN]))
                .setChromosomeBand(values[CHROMOSOME_BAND_COLUMN])
                .setEntrezId(entrezIds)
                .setGene(gene)
                .setGeneID(values[GENE_ID_COLUMN])
                .setGeneStart(Long.valueOf(values[GENE_START_COLUMN]))
                .setGeneEnd(Long.valueOf(values[GENE_END_COLUMN]))
                .setCodingStart(codingStart)
                .setCodingEnd(codingEnd);
    }
}
