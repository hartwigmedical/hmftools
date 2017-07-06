package com.hartwig.hmftools.common.zipper;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@Deprecated
public class GenomeZipper<R extends GenomeRegion> {

    public GenomeZipper(@NotNull final List<R> regions, @NotNull final GenomeZipperRegionHandler<R> regionHandler) {
        this.regions = regions;
        this.regionHandler = regionHandler;
    }

    @NotNull
    private final List<R> regions;
    @NotNull
    private final GenomeZipperRegionHandler<R> regionHandler;
    private final List<Consumer<GenomeRegion>> regionConsumers = Lists.newArrayList();

    public <P extends GenomePosition> void addPositions(@NotNull final List<P> positions, @NotNull final Consumer<P> handler) {
        regionConsumers.add(new PositionConsumer<>(GenomePositionSelectorFactory.create(positions), handler));
    }

    public void zip() {
        String chromosome = "";

        for (final R region : regions) {

            if (!region.chromosome().equals(chromosome)) {
                chromosome = region.chromosome();
                regionHandler.chromosome(chromosome);
            }

            regionHandler.enter(region);
            for (Consumer<GenomeRegion> regionConsumer : regionConsumers) {
                regionConsumer.accept(region);
            }
            regionHandler.exit(region);
        }
    }

    class PositionConsumer<P extends GenomePosition> implements Consumer<GenomeRegion> {
        private final GenomePositionSelector<P> selector;
        private final Consumer<P> handler;

        PositionConsumer(final GenomePositionSelector<P> selector, final Consumer<P> handler) {
            this.selector = selector;
            this.handler = handler;
        }

        @Override
        public void accept(final GenomeRegion genomeRegion) {
            selector.select(genomeRegion, handler);
        }
    }
}