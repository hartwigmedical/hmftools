package com.hartwig.hmftools.common.purple.region;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum ConsolidatedRegionWriter {
    ;

    public static void writeRegions(@NotNull final String filePath, @NotNull Collection<ConsolidatedRegion> purity)
            throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(ConsolidatedRegionWriter::transform).forEach(lines::add);

        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return "chromosome" + '\t' + "start" + '\t' + "end" + '\t' + "copyNumber" + '\t' + "bafCount" + '\t'
                + "observedBAF" + '\t' + "actualBAF";
    }

    @NotNull
    private static String transform(@NotNull final ConsolidatedRegion region) {
        return region.chromosome() + '\t' + region.start() + '\t' + region.end() + '\t'
                + region.averageTumorCopyNumber() + '\t' + region.bafCount() + '\t' + region.averageObservedBAF()
                + '\t' + region.averageActualBAF();
    }
}
