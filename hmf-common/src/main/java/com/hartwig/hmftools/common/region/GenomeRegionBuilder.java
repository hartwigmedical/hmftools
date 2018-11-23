package com.hartwig.hmftools.common.region;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class GenomeRegionBuilder {

    private final String chromosome;
    private final int minGap;
    @NotNull
    private final List<GenomeRegion> regions;

    public GenomeRegionBuilder(@NotNull final String chromosome) {
        this(chromosome, 1);
    }

    public GenomeRegionBuilder(@NotNull final String chromosome, final int minGap) {
        this.chromosome = chromosome;
        this.minGap = minGap;
        this.regions = Lists.newArrayList();
    }

    @NotNull
    public List<GenomeRegion> build() {
        return regions;
    }

    public void addPosition(final long position) {
        GenomeRegion prev = null;

        for (int i = 0; i < regions.size(); i++) {
            GenomeRegion current = regions.get(i);
            if (position >= current.start() && position <= current.end()) {
                return;
            }

            if (position < current.start()) {
                // Attach to previous?
                if (prev != null && prev.end() + minGap >= position) {

                    if (position + minGap >= current.start()) {
                        // Join previous and current
                        prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), current.end());
                        regions.set(i - 1, prev);
                        regions.remove(i);
                        return;
                    } else {
                        // Add to previous
                        prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), position);
                        regions.set(i - 1, prev);
                        return;
                    }
                }

                // Attach to current
                if (position + minGap >= current.start()) {
                    current = GenomeRegionFactory.create(current.chromosome(), position, current.end());
                    regions.set(i, current);
                    return;
                }

                // Attach between
                regions.add(i, GenomeRegionFactory.create(chromosome, position, position));
                return;
            }

            prev = current;
        }

        if (prev != null && prev.end() + minGap >= position) {
            prev = GenomeRegionFactory.create(prev.chromosome(), prev.start(), position);
            regions.set(regions.size() - 1, prev);
        } else {
            regions.add(GenomeRegionFactory.create(chromosome, position, position));
        }

    }

}
