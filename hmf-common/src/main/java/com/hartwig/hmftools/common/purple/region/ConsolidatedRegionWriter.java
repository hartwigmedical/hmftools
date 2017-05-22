package com.hartwig.hmftools.common.purple.region;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum ConsolidatedRegionWriter {
    ;

    public static void writeRegions(@NotNull final String filePath, @NotNull Collection<ConsolidatedRegion> purity) throws IOException {

        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(ConsolidatedRegionWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header() {
        return new StringBuilder()
                .append("chromosome").append('\t')
                .append("start").append('\t')
                .append("end").append('\t')
                .append("copyNumber").append('\t')
                .append("bafCount")
                .append("observedBAF")
                .toString();
    }

    private static String transform(ConsolidatedRegion region) {
        return new StringBuilder()
                .append(region.chromosome()).append('\t')
                .append(region.start()).append('\t')
                .append(region.end()).append('\t')
                .append(region.averageTumorCopyNumber()).append('\t')
                .append(region.bafCount())
                .append(region.averageObservedBAF())
                .toString();
    }

}
