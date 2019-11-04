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
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public class HmfGenePanelBed {

    private static final String DELIMITER = "\t";
    static final int EXTRA_BASES = 2;

    public static void write37File(@NotNull final String filename) throws IOException {
        writeBedFile(filename, createRegions(HmfGenePanelSupplier.allGeneList37()));
    }

    public static void write38File(@NotNull final String filename) throws IOException {
        writeBedFile(filename, createRegions(HmfGenePanelSupplier.allGeneList38()));
    }

    private static void writeBedFile(@NotNull final String filename, @NotNull final List<GenomeRegion> regions) throws IOException {
        List<String> strings = regions.stream().map(HmfGenePanelBed::asBed).collect(Collectors.toList());
        Files.write(new File(filename).toPath(), strings);
    }

    @NotNull
    static List<GenomeRegion> createRegions(@NotNull final List<HmfTranscriptRegion> regions) {
        final Map<String, GenomeRegions> regionsMap = Maps.newHashMap();

        for (HmfTranscriptRegion region : regions) {
            final GenomeRegions regionBuilder = regionsMap.computeIfAbsent(region.chromosome(), GenomeRegions::new);
            regionBuilder.addRegion(region.start() - EXTRA_BASES, region.end() + EXTRA_BASES);
        }

        final List<GenomeRegion> result = Lists.newArrayList();
        for (GenomeRegions genomeRegions : regionsMap.values()) {
            result.addAll(genomeRegions.build());
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    private static String asBed(@NotNull final GenomeRegion region) {
        return new StringJoiner(DELIMITER).add(String.valueOf(region.chromosome()))
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }

}
