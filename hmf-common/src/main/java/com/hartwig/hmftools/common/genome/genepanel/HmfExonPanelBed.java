package com.hartwig.hmftools.common.genome.genepanel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

// This class is only used by tests but we like to keep for future reference
final class HmfExonPanelBed {

    private static final String DELIMITER = "\t";

    @SuppressWarnings("WeakerAccess")
    public static void write19File(@NotNull final String filename) throws IOException {
        DriverGenePanel genePanel = new DriverGenePanelFactory().create();
        writeBedFile(filename, Strings.EMPTY, createRegions(genePanel, HmfGenePanelSupplier.allGeneList37()));
    }

    @SuppressWarnings("WeakerAccess")
    public static void write38File(@NotNull final String filename) throws IOException {
        DriverGenePanel genePanel = new DriverGenePanelFactory().create();
        writeBedFile(filename, "chr", createRegions(genePanel, HmfGenePanelSupplier.allGeneList38()));
    }

    private static void writeBedFile(@NotNull final String filename, @NotNull final String prefix,
            @NotNull final List<GenomeRegion> regions) throws IOException {
        List<String> strings = regions.stream().map(x -> asBed(prefix, x)).collect(Collectors.toList());
        Files.write(new File(filename).toPath(), strings);
    }

    @NotNull
    static List<GenomeRegion> createRegions(@NotNull final DriverGenePanel genePanel, @NotNull final List<HmfTranscriptRegion> regions) {
        final Set<String> actionableGenes = Sets.newHashSet();
        actionableGenes.addAll(genePanel.oncoGenes());
        actionableGenes.addAll(genePanel.tsGenes());

        final Map<String, GenomeRegions> regionsMap = Maps.newHashMap();

        for (HmfTranscriptRegion transcript : regions) {
            final GenomeRegions regionBuilder = regionsMap.computeIfAbsent(transcript.chromosome(), GenomeRegions::new);

            boolean forward = transcript.strand() == Strand.FORWARD;
            for (int i = 0; i < transcript.exome().size(); i++) {
                if (transcript.codingStart() == 0 || !actionableGenes.contains(transcript.gene())) {
                    continue;
                }

                // Splice sites (+1,+2,+5, -2,-1)
                final HmfExonRegion exon = transcript.exome().get(i);
                if (i != 0) {
                    regionBuilder.addRegion(exon.start() - 2, exon.start() - 1);
                    if (!forward) {
                        regionBuilder.addPosition(exon.start() - 5);
                    }
                }
                if (i != transcript.exome().size() - 1) {
                    regionBuilder.addRegion(exon.end() + 1, exon.end() + 2);
                    if (forward) {
                        regionBuilder.addPosition(exon.end() + 5);
                    }
                }

                if (transcript.codingStart() < exon.end() && transcript.codingEnd() > exon.start()) {
                    regionBuilder.addRegion(Math.max(transcript.codingStart(), exon.start()), Math.min(transcript.codingEnd(), exon.end()));
                }
            }
        }

        final List<GenomeRegion> result = Lists.newArrayList();
        for (GenomeRegions genomeRegions : regionsMap.values()) {
            result.addAll(genomeRegions.build());
        }

        Collections.sort(result);

        return result;
    }

    @NotNull
    private static String asBed(@NotNull final String prefix, @NotNull final GenomeRegion region) {
        return new StringJoiner(DELIMITER).add(prefix + region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }
}
