package com.hartwig.hmftools.svanalysis.visualisation;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

public class Span {

    @NotNull
    public static List<GenomeRegion> span(@NotNull final List<Segment> segments, @NotNull final List<Link> links) {

        final List<GenomePosition> positions = Lists.newArrayList();
        positions.addAll(Segments.allPositions(segments));
        positions.addAll(Links.allPositions(links));

        return span(positions);
    }

    @NotNull
    public static List<GenomeRegion> span(@NotNull final List<GenomePosition> positions) {
        final List<GenomeRegion> result = Lists.newArrayList();

        final List<String> chromosomes = positions.stream().map(GenomePosition::chromosome).distinct().collect(Collectors.toList());
        for (final String chromosome : chromosomes) {
            long min =
                    positions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomePosition::position).min().orElse(0);
            long max =
                    positions.stream().filter(x -> x.chromosome().equals(chromosome)).mapToLong(GenomePosition::position).max().orElse(0);

            result.add(GenomeRegionFactory.create(chromosome, min, max));
        }

        Collections.sort(result);
        return result;
    }

}
