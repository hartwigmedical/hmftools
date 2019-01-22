package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class CopyNumberAlterations {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";


    @NotNull
    public static List<CopyNumberAlteration> copyNumberInTracks(@NotNull final List<CopyNumberAlteration> raw, @NotNull final List<Track> tracks) {
        final List<CopyNumberAlteration> result = Lists.newArrayList();

        for (CopyNumberAlteration copyNumberAlteration : raw) {
            final String contig = copyNumberAlteration.chromosome();
            final List<Track> chromosomeTracks = tracks.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if (!chromosomeTracks.isEmpty()) {
                long minPosition = chromosomeTracks.stream().mapToLong(GenomeRegion::start).min().orElse(0);
                long maxPosition = chromosomeTracks.stream().mapToLong(GenomeRegion::end).min().orElse(0);
                if (copyNumberAlteration.end() >= minPosition && copyNumberAlteration.start() <= maxPosition) {
                    result.add(ImmutableCopyNumberAlteration.builder()
                            .from(copyNumberAlteration)
                            .start(Math.max(minPosition, copyNumberAlteration.start()))
                            .end(Math.min(maxPosition, copyNumberAlteration.end()))
                            .build());
                }
            }

        }

        return result;
    }

    @NotNull
    public static List<CopyNumberAlteration> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<CopyNumberAlteration> fromLines(@NotNull List<String> lines) {
        final List<CopyNumberAlteration> results = Lists.newArrayList();
        for (final String line : lines) {
            if (!line.startsWith(COMMENT)) {
                results.add(fromString(line));
            }
        }

        return results;
    }

    @NotNull
    private static CopyNumberAlteration fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableCopyNumberAlteration.builder()
                .chromosome(values[0])
                .start(Long.valueOf(values[1]))
                .end(Long.valueOf(values[2]))
                .copyNumber(Double.valueOf(values[3]))
                .build();
    }

}
