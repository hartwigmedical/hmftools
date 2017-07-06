package com.hartwig.hmftools.common.zipper;

import java.util.List;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Deprecated
public class SimpleGenomeZipper<R extends GenomeRegion, P extends GenomePosition>
        implements GenomeZipperRegionHandler<R> {

    public static <R extends GenomeRegion, P extends GenomePosition> void zip(@NotNull List<R> regions, @NotNull List<P> positions,
            @NotNull SimpleGenomeZipperInRegionPositionsHandler<R, P> handler) {
        zipInner(regions, positions, handler::handle);
    }

    private static <R extends GenomeRegion, P extends GenomePosition> void zipInner( @NotNull List<R> regions, @NotNull List<P> positions, @NotNull BiConsumer<R, P> handler) {
        SimpleGenomeZipper<R, P> simpleZipper = new SimpleGenomeZipper<>(handler);
        GenomeZipper<R> genomeZipper = new GenomeZipper<>(regions, simpleZipper);
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
    public void chromosome(@NotNull final String chromosome) {

    }

    @Override
    public void enter(@NotNull final R region) {
        this.region = region;
    }

    @Override
    public void exit(@NotNull final R region) {
        this.region = null;
    }

    private void handlePosition(P position) {
        handler.accept(region, position);
    }
}
