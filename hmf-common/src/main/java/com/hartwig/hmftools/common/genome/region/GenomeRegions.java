package com.hartwig.hmftools.common.genome.region;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class GenomeRegions {

    @NotNull
    public static GenomeRegion create(@NotNull final String chromosome, final int start, final int end) {
        return ImmutableGenomeRegionImpl.builder().chromosome(chromosome).start(start).end(end).build();
    }

    private final String chromosome;
    private final int minGap;
    @NotNull
    private final List<GenomeRegion> regions;

    public GenomeRegions(@NotNull final String chromosome) {
        this(chromosome, 1);
    }

    public GenomeRegions(@NotNull final String chromosome, final int minGap) {
        this.chromosome = chromosome;
        this.minGap = minGap;
        this.regions = Lists.newArrayList();
    }

    @NotNull
    public List<GenomeRegion> build() {
        return regions;
    }

    public void addRegion(final int start, final int end) {
        addPosition(start);
        int i = indexOf(start);

        final GenomeRegion current = regions.get(i);
        int extendedEnd = Math.max(current.end(), end);

        while (i + 1 < regions.size()) {
            final GenomeRegion next = regions.get(i + 1);
            if (next.start() <= extendedEnd + minGap) {
                regions.remove(i + 1);
                extendedEnd = Math.max(extendedEnd, next.end());
            } else {
                break;
            }
        }

        final GenomeRegion extendedRegion = GenomeRegions.create(current.chromosome(), current.start(), extendedEnd);
        regions.set(i, extendedRegion);
    }

    public void addPosition(final int position) {
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
                        prev = create(prev.chromosome(), prev.start(), current.end());
                        regions.set(i - 1, prev);
                        regions.remove(i);
                        return;
                    } else {
                        // Add to previous
                        prev = create(prev.chromosome(), prev.start(), position);
                        regions.set(i - 1, prev);
                        return;
                    }
                }

                // Attach to current
                if (position + minGap >= current.start()) {
                    current = create(current.chromosome(), position, current.end());
                    regions.set(i, current);
                    return;
                }

                // Attach between
                regions.add(i, create(chromosome, position, position));
                return;
            }

            prev = current;
        }

        if (prev != null && prev.end() + minGap >= position) {
            prev = GenomeRegions.create(prev.chromosome(), prev.start(), position);
            regions.set(regions.size() - 1, prev);
        } else {
            regions.add(create(chromosome, position, position));
        }

    }

    private int indexOf(int position) {
        for (int i = regions.size() - 1; i >= 0; i--) {
            GenomeRegion region = regions.get(i);
            if (position >= region.start() && position <= region.end()) {
                return i;
            }
        }

        return -1;
    }

    /**
     * Create genome regions from a sorted list of genome positions. This version is much faster then using
     * GenomeRegionsBuilder and add positions one by one.
     *
     * @param minGap                  minimum gap between the returned genome regions.
     * @param sortedGenomePositions   input list must be sorted.
     * @return genome regions         list of genome region that contains all input genome positions.
     */
    public static List<GenomeRegion> fromSortedGenomePositions(int minGap, final List<? extends GenomePosition> sortedGenomePositions)
    {
        List<GenomeRegion> genomeRegions = new ArrayList<>();

        String chromosome = null;
        int start = -1;
        int end = -1;

        for (GenomePosition position : sortedGenomePositions)
        {
            if (chromosome == null || !chromosome.equals(position.chromosome()) || end + minGap < position.position())
            {
                if (chromosome != null)
                {
                    // we finalise previous one
                    genomeRegions.add(GenomeRegions.create(chromosome, start, end));
                }

                chromosome = position.chromosome();
                start = position.position();
                end = position.position();
                continue;
            }

            if (end > position.position())
            {
                // this means the input genome positions are not sorted
                throw new IllegalArgumentException("Genome position going backwards, input might not be sorted");
            }

            // add to existing region
            end = position.position();
        }

        if (chromosome != null)
        {
            // we finalise last one
            genomeRegions.add(GenomeRegions.create(chromosome, start, end));
        }

        return genomeRegions;
    }
}
