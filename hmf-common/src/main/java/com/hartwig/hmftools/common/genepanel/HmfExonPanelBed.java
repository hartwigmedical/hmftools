package com.hartwig.hmftools.common.genepanel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.common.region.HmfExonRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class HmfExonPanelBed {

    private static final String DELIMITER = "\t";
    static final int EXTRA_BASES = 2;

    public static void write19File(@NotNull final String filename) throws IOException {
        writeBedFile(filename, Strings.EMPTY, createRegions(HmfGenePanelSupplier.allGeneList37()));
    }

    public static void write38File(@NotNull final String filename) throws IOException {
        writeBedFile(filename, "chr", createRegions(HmfGenePanelSupplier.allGeneList38()));
    }

    private static void writeBedFile(@NotNull final String filename, @NotNull final String prefix,
            @NotNull final List<GenomeRegion> regions) throws IOException {
        List<String> strings = regions.stream().map(x -> asBed(prefix, x)).collect(Collectors.toList());
        Files.write(new File(filename).toPath(), strings);
    }

    @NotNull
    static List<GenomeRegion> createRegions(@NotNull final List<HmfTranscriptRegion> regions) {
        final Map<String, GenomeRegions> regionsMap = Maps.newHashMap();

        for (HmfTranscriptRegion transcript : regions) {
            for (HmfExonRegion exon : transcript.exome()) {
                final GenomeRegions regionBuilder = regionsMap.computeIfAbsent(exon.chromosome(), GenomeRegions::new);
                regionBuilder.addRegion(exon.start() - EXTRA_BASES, exon.end() + EXTRA_BASES);
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
        return new StringJoiner(DELIMITER).add(prefix + String.valueOf(region.chromosome()))
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }

}
