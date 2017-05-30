package com.hartwig.hmftools.common.zipper;

import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class GenomeZipper<R extends GenomeRegion> {

    public GenomeZipper(@NotNull final List<R> regions, @NotNull final GenomeZipperRegionHandler<R> regionHandler) {
        this.regions = regions;
        this.regionHandler = regionHandler;
    }

    @NotNull
    private final List<R> regions;
    @NotNull
    private final GenomeZipperRegionHandler<R> regionHandler;
    private final List<InnerPosition<?>> innerPositions = Lists.newArrayList();

    public <P extends GenomePosition> void addPositions(@NotNull List<P> positions, @NotNull Consumer<P> handler) {
        addPositions(positions, (x, y) -> handler.accept(y));
    }

    public <P extends GenomePosition> void addPositions(@NotNull List<P> positions,
            @NotNull BiConsumer<R, P> handler) {
        innerPositions.add(new InnerPosition<>(positions, handler));
    }

    public void run() {
        for (final R region : regions) {
            regionHandler.enter(region);
            for (final InnerPosition<?> innerPosition : innerPositions) {
                innerPosition.apply(region);
            }
            regionHandler.exit(region);
        }
    }

    private class InnerPosition<P extends GenomePosition> {

        @NotNull
        private final List<P> positions;
        @NotNull
        private final BiConsumer<R, P> handler;

        private int index;

        private InnerPosition(@NotNull final List<P> positions, @NotNull final BiConsumer<R, P> handler) {
            this.positions = positions;
            this.handler = handler;
        }

        private void apply(R region) {
            while (index < positions.size()) {
                P position = positions.get(index);
                int compare = compare(position, region);
                if (compare > 0) {
                    break;
                } else if (compare == 0) {
                    handler.accept(region, position);
                }

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