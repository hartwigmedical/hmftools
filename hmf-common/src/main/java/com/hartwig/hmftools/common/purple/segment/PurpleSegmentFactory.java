package com.hartwig.hmftools.common.purple.segment;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.zipper.GenomeZipper;
import com.hartwig.hmftools.common.zipper.GenomeZipperRegionHandler;

import org.jetbrains.annotations.NotNull;

public class PurpleSegmentFactory implements GenomeZipperRegionHandler<GenomeRegion> {

    public static List<PurpleSegment> createSegmentsOuter(List<GenomeRegion> regions, List<StructuralVariant> variants) {
        return createSegments(regions, StructuralVariantSegments.create(variants));
    }

    @VisibleForTesting
    static List<PurpleSegment> createSegments(List<GenomeRegion> regions, List<StructuralVariantSegment> variants) {
        final PurpleSegmentFactory factory = new PurpleSegmentFactory();
        final GenomeZipper<GenomeRegion> zipper = new GenomeZipper<>(false, regions, factory);
        zipper.addPositions(variants, factory::variant);
        zipper.zip();
        return factory.segments();
    }

    @VisibleForTesting
    static final long MIN_BASES = 2001;
    private final List<PurpleSegment> segments;

    private String chromosome = null;
    private long start;
    private long end;

    private PurpleSegmentFactory() {
        segments = Lists.newArrayList();
    }

    @NotNull
    private PurpleSegmentSource source = PurpleSegmentSource.FREEC;

    @NotNull
    private List<PurpleSegment> segments() {
        return segments;
    }

    @Override
    public void enter(final GenomeRegion region) {
        chromosome = region.chromosome();
        start = region.start();
        end = region.end();
        source = PurpleSegmentSource.FREEC;
    }

    @Override
    public void exit(final GenomeRegion region) {
        insert(region.end());
    }

    public void variant(final StructuralVariantSegment variant) {

        if (variant.position() - start < MIN_BASES) {
            // do nothing
        } else if (end - variant.position() < MIN_BASES) {
            // do nothing
        } else {
            insert(variant.position() - 1);
            start = variant.position();
            source = PurpleSegmentSource.fromVariant(variant.type());
        }
    }

    private void insert(long end) {
        segments.add(ImmutablePurpleSegment.builder().chromosome(chromosome).start(start).end(end).source(source).build());
    }
}
