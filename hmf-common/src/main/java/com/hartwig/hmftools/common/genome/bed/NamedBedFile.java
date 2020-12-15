package com.hartwig.hmftools.common.genome.bed;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.BEDFileLoader;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;

public final class NamedBedFile {

    private static final Logger LOGGER = LogManager.getLogger(BEDFileLoader.class);
    private static final String DELIMITER = "\t";

    public static void writeBedFile(@NotNull final String filename, @NotNull final List<NamedBed> regions) throws IOException {
        List<String> strings = regions.stream().map(NamedBedFile::asBed).collect(Collectors.toList());
        Files.write(new File(filename).toPath(), strings);
    }

    @NotNull
    public static List<NamedBed> readBedFile(@NotNull String bedFile) throws IOException {
        List<NamedBed> result = Lists.newArrayList();
        NamedBed prevRegion = null;
        try (final AbstractFeatureReader<BEDFeature, LineIterator> reader = getFeatureReader(bedFile, new BEDCodec(), false)) {
            for (final BEDFeature bedFeature : reader.iterator()) {
                final NamedBed namedBed = fromBedFeature(bedFeature);
                if (namedBed.end() < namedBed.start()) {
                    LOGGER.warn("Invalid genome region found in chromosome {}: start={}, end={}",
                            namedBed.chromosome(),
                            namedBed.start(),
                            namedBed.end());
                } else {
                    if (prevRegion != null && namedBed.chromosome().equals(prevRegion.chromosome())
                            && prevRegion.end() >= namedBed.start()) {
                        LOGGER.warn("BED file is not sorted, please fix! Current={}, Previous={}", namedBed, prevRegion);
                    } else {
                        result.add(namedBed);
                        prevRegion = namedBed;
                    }
                }
            }
        }

        return result;
    }

    @NotNull
    static NamedBed fromBedFeature(@NotNull BEDFeature feature) {
        String name = feature.getName();
        return ImmutableNamedBed.builder()
                .chromosome(feature.getContig())
                .start(feature.getStart())
                .end(feature.getEnd())
                .name(name == null ? Strings.EMPTY : name)
                .build();
    }

    @NotNull
    private static String asBed(@NotNull final NamedBed region) {
        return new StringJoiner(DELIMITER).add(region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .add(region.name())
                .toString();
    }
}
