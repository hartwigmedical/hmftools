package com.hartwig.hmftools.common.zipper;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class GenomeZipper<R extends GenomeRegion> {

    private final boolean includePositionsOutsideRegions;

    public GenomeZipper(final boolean includePositionsOutsideRegions, @NotNull final List<R> regions, @NotNull final GenomeZipperRegionHandler<R> regionHandler) {
        this.includePositionsOutsideRegions = includePositionsOutsideRegions;
        this.regions = regions;
        this.regionHandler = regionHandler;
    }

    @NotNull
    private final List<R> regions;
    @NotNull
    private final GenomeZipperRegionHandler<R> regionHandler;
    private final List<InnerPosition<?>> innerPositions = Lists.newArrayList();

    public <P extends GenomePosition> void addPositions(@NotNull final List<P> positions,
            @NotNull final Consumer<P> handler) {
        innerPositions.add(new InnerPosition<>(positions, handler));
    }

    public void zip() {
        String chromosome = "";

        for (final R region : regions) {

            if (!region.chromosome().equals(chromosome)) {
                chromosome = region.chromosome();
                regionHandler.chromosome(chromosome);
            }

            if (includePositionsOutsideRegions) {
                for (final InnerPosition<?> innerPosition : innerPositions) {
                    innerPosition.preRegion(region);
                }
            }

            regionHandler.enter(region);
            for (final InnerPosition<?> innerPosition : innerPositions) {
                innerPosition.inRegion(region);
            }
            regionHandler.exit(region);
        }

        if (includePositionsOutsideRegions) {
            for (final InnerPosition<?> innerPosition : innerPositions) {
                innerPosition.remainder();
            }
        }
    }

    private class InnerPosition<P extends GenomePosition> {

        @NotNull
        private final List<P> positions;
        @NotNull
        private final Consumer<P> handler;

        private int index;

        private InnerPosition(@NotNull final List<P> positions, @NotNull final Consumer<P> handler) {
            this.positions = positions;
            this.handler = handler;
        }

        private void preRegion(@NotNull final R region) {
            while (index < positions.size()) {
                final P position = positions.get(index);
                final int compare = compare(position, region);
                if (compare >= 0) {
                    break;
                } else{
                    handler.accept(position);
                }

                index++;
            }
        }

        private void inRegion(@NotNull final R region) {
            while (index < positions.size()) {
                final P position = positions.get(index);
                final int compare = compare(position, region);
                if (compare > 0) {
                    break;
                } else if (compare == 0) {
                    handler.accept(position);
                }

                index++;
            }
        }

        private void remainder() {
            while (index < positions.size()) {
                final P position = positions.get(index);
                handler.accept(position);
                index++;
            }
        }
    }

    private static int compare(@NotNull final GenomePosition position, @NotNull final GenomeRegion region) {
        int positionChromosome = Chromosomes.asInt(position.chromosome());
        int regionChromosome = Chromosomes.asInt(region.chromosome());
        if (positionChromosome < regionChromosome) {
            return -1;
        }
        if (positionChromosome > regionChromosome) {
            return 1;
        }

        if (position.position() < region.start()) {
            return -1;
        }

        if (position.position() > region.end()) {
            return 1;
        }

        return 0;
    }
}