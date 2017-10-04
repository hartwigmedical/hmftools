package com.hartwig.hmftools.purple.segment;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.pcf.ImmutablePCFRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class SegmentMerge implements RegionZipperHandler<GenomeRegion, GenomeRegion> {

    static List<GenomeRegion> merge(List<GenomeRegion> normal, List<GenomeRegion> tumor) {
        SegmentMerge handler = new SegmentMerge();
        RegionZipper.zip(normal, tumor, handler);
        handler.enterChromosome("finalise");
        return handler.regions;
    }

    private String chromosome;
    private GenomeRegion previous;
    private List<GenomeRegion> regions = Lists.newArrayList();

    @Override
    public void enterChromosome(@NotNull final String chromosome) {
        if (previous != null) {
            addRegion(previous.start(), previous.end());
        }
        previous = null;
        this.chromosome = chromosome;
    }

    @Override
    public void primary(@NotNull final GenomeRegion region) {
        previous = merge(previous, region);
    }

    @Override
    public void secondary(@NotNull final GenomeRegion region) {
        primary(region);
    }

    @Nullable
    private GenomeRegion merge(@Nullable GenomeRegion first, @NotNull GenomeRegion second) {
        if (first == null) {
            return second;
        }

        if (first.end() < second.start()) {
            addRegion(first.start(), first.end());
            return second;
        }

        if (first.start() == second.start()) {
            if (second.end() < first.end()) {
                addRegion(first.start(), second.end());
                return create(second.end() + 1, first.end());
            } else if (second.end() == first.end()) {
                addRegion(first.start(), first.end());
                return null;
            } else {
                addRegion(first.start(), first.end());
                return create(first.end() + 1, second.end());
            }
        }

        addRegion(first.start(), second.start() - 1);
        return merge(create(second.start(), first.end()), second);
    }

    private void addRegion(long start, long end) {
        regions.add(create(start, end));
    }

    private GenomeRegion create(long start, long end) {
        return ImmutablePCFRegion.builder().chromosome(chromosome).start(start).end(end).build();
    }
}
