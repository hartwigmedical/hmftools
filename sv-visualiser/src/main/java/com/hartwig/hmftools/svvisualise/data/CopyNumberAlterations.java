package com.hartwig.hmftools.svvisualise.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class CopyNumberAlterations {



    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<CopyNumberAlteration> copyNumbers(long copyNumberDistance, @NotNull final List<CopyNumberAlteration> alterations,
            @NotNull final List<GenomeRegion> span) {
        final List<CopyNumberAlteration> result = Lists.newArrayList();

        for (int i = 0; i < alterations.size(); i++) {
            CopyNumberAlteration alteration = alterations.get(i);
            final String contig = alteration.chromosome();
            final List<GenomeRegion> chromosomeSegments =
                    span.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if (!chromosomeSegments.isEmpty()) {
                long minTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::start).min().orElse(0) - copyNumberDistance;
                long maxTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::end).max().orElse(0) + copyNumberDistance;
                if (alteration.end() >= minTrackPosition && alteration.start() <= maxTrackPosition) {

                    boolean isStartDecreasing = i > 0 && lessThan(alteration, alterations.get(i - 1));
                    long startPosition = isStartDecreasing ? alteration.start() - 1 : alteration.start();

                    boolean isEndIncreasing = i < alterations.size() - 1 && lessThan(alteration, alterations.get(i + 1));
                    long endPosition = isEndIncreasing ? alteration.end() + 1 : alteration.end();

                    result.add(ImmutableCopyNumberAlteration.builder()
                            .from(alteration)
                            .start(Math.max(minTrackPosition, startPosition))
                            .end(Math.min(maxTrackPosition, endPosition))
                            .build());
                }
            }

        }

        return result;
    }

    private static boolean lessThan(@NotNull final CopyNumberAlteration first, @NotNull final CopyNumberAlteration second) {
        return first.chromosome().equals(second.chromosome()) && Doubles.lessThan(first.copyNumber(), second.copyNumber());
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
                .sampleId(values[0])
                .chromosome(values[1])
                .start(Long.valueOf(values[2]))
                .end(Long.valueOf(values[3]))
                .copyNumber(Double.valueOf(values[4]))
                .baf(Double.valueOf(values[5]))
                .build();
    }

}
