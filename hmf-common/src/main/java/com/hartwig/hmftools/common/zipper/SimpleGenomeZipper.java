package com.hartwig.hmftools.common.zipper;

import java.util.List;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SimpleGenomeZipper<R extends GenomeRegion, P extends GenomePosition>
        implements GenomeZipperRegionHandler<R> {

    public static <R extends GenomeRegion, P extends GenomePosition> void zip(@NotNull List<R> regions, @NotNull List<P> positions,
            @NotNull SimpleGenomeZipperAllPositionsHandler<R, P> handler) {
        zipInner(true, regions, positions, handler::handle);
    }

    public static <R extends GenomeRegion, P extends GenomePosition> void zip(@NotNull List<R> regions, @NotNull List<P> positions,
            @NotNull SimpleGenomeZipperInRegionPositionsHandler<R, P> handler) {
        zipInner(false, regions, positions, handler::handle);
    }

    private static <R extends GenomeRegion, P extends GenomePosition> void zipInner(
            boolean includePositionsOutsideRegion, @NotNull List<R> regions, @NotNull List<P> positions, @NotNull BiConsumer<R, P> handler) {
        SimpleGenomeZipper<R, P> simpleZipper = new SimpleGenomeZipper<>(handler);
        GenomeZipper<R> genomeZipper = new GenomeZipper<>(includePositionsOutsideRegion, regions, simpleZipper);
        genomeZipper.addPositions(positions, simpleZipper::handlePosition);
        genomeZipper.zip();
    }

    private final BiConsumer<R, P> handler;

    @Nullable
    private R region;

    private SimpleGenomeZipper(final BiConsumer<R, P> handler) {
        this.handler = handler;
    }

    @Override
    public void enter(final R region) {
        this.region = region;
    }

    @Override
    public void exit(final R region) {
        this.region = null;
    }

    private void handlePosition(P position) {
        handler.accept(region, position);
    }
}
