package com.hartwig.hmftools.common.purple.segment;

import static com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport.NONE;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.zipper.GenomeZipper;
import com.hartwig.hmftools.common.zipper.GenomeZipperRegionHandler;

import org.jetbrains.annotations.NotNull;

public class PurpleSegmentFactory implements GenomeZipperRegionHandler<GenomeRegion> {

    public static List<PurpleSegment> createSegments(List<GenomeRegion> regions, List<StructuralVariant> variants) {
        return createSegmentsInner(regions, StructuralVariantPositions.create(variants));
    }

    public static List<PurpleSegment> createSegments(List<GenomeRegion> regions) {
        return regions.stream()
                .map(x -> ImmutablePurpleSegment.builder().from(x).ratioSupport(true).structuralVariantSupport(NONE).build())
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static List<PurpleSegment> createSegmentsInner(List<GenomeRegion> regions, List<StructuralVariantPosition> variants) {
        final PurpleSegmentFactory factory = new PurpleSegmentFactory();
        final GenomeZipper<GenomeRegion> zipper = new GenomeZipper<>(false, regions, factory);
        zipper.addPositions(variants, factory::variant);
        zipper.zip();
        return factory.segments();
    }

    @VisibleForTesting
    static final long MIN_BASES = 1001;
    private final List<PurpleSegment> segments;

    private String chromosome = null;
    private long start;
    private long end;

    private PurpleSegmentFactory() {
        segments = Lists.newArrayList();
    }

    @NotNull
    private StructuralVariantSupport svStart = NONE;
    @NotNull
    private StructuralVariantSupport svEnd = NONE;
    private boolean ratioSupport;

    @NotNull
    private List<PurpleSegment> segments() {
        return segments;
    }

    @Override
    public void chromosome(@NotNull final String chromosome) {
        this.chromosome = chromosome;
        svStart = NONE;
        svEnd = NONE;
    }

    @Override
    public void enter(@NotNull final GenomeRegion region) {
        start = region.start();
        end = region.end();
        ratioSupport = true;
    }

    @Override
    public void exit(@NotNull final GenomeRegion region) {
        insert(region.end());
        if (!svEnd.equals(NONE)) {
            svStart = svEnd;
        } else {
            svStart = NONE;
        }
        svEnd = NONE;
    }

    public void variant(final StructuralVariantPosition variant) {

        if (variant.position() - start < MIN_BASES) {
            svStart = svSupport(svStart, variant);
        } else if (end - variant.position() < MIN_BASES) {
            svEnd = svSupport(svEnd, variant);
        } else {
            insert(variant.position() - 1);
            ratioSupport = false;
            start = variant.position();
            svStart = StructuralVariantSupport.fromVariant(variant.type());
        }
    }

    private StructuralVariantSupport svSupport(StructuralVariantSupport current, StructuralVariantPosition variant) {
        return current.equals(NONE) ? StructuralVariantSupport.fromVariant(variant.type()) : MULTIPLE;
    }

    private void insert(long end) {
        segments.add(ImmutablePurpleSegment.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .ratioSupport(ratioSupport)
                .structuralVariantSupport(svStart)
                .build());
    }
}
