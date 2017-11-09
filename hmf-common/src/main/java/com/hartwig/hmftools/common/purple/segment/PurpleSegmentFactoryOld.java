package com.hartwig.hmftools.common.purple.segment;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

@Deprecated
public class PurpleSegmentFactoryOld {

    public static List<PurpleSegment> createSegments(List<GenomeRegion> regions, List<StructuralVariant> variants) {
        return createSegmentsInner(regions, StructuralVariantPositionFactory.createMinBaseFiltered(variants));
    }

    @VisibleForTesting
    static List<PurpleSegment> createSegmentsInner(List<GenomeRegion> regions, List<StructuralVariantPosition> variants) {
        return new PurpleSegmentFactoryOld(regions, variants).segments();
    }

    @VisibleForTesting
    static final long MIN_BASES = 1001;
    private final List<PurpleSegment> segments;

    private String chromosome = null;
    private long start;
    private long end;

    private PurpleSegmentFactoryOld(List<GenomeRegion> regions, List<StructuralVariantPosition> variants) {
        segments = Lists.newArrayList();

        final GenomePositionSelector<StructuralVariantPosition> variantSelector = GenomePositionSelectorFactory.create(variants);
        String chromosome = "";
        for (final GenomeRegion region : regions) {

            if (!chromosome.equals(region.chromosome())) {
                chromosome = region.chromosome();
                chromosome(chromosome);
            }

            enter(region);
            variantSelector.select(region, this::variant);
            exit(region);
        }
    }

    @NotNull
    private SegmentSupport svStart = NONE;
    @NotNull
    private SegmentSupport svEnd = NONE;
    private boolean ratioSupport;

    @NotNull
    private List<PurpleSegment> segments() {
        return segments;
    }

    private void chromosome(@NotNull final String chromosome) {
        this.chromosome = chromosome;
        svStart = NONE;
        svEnd = NONE;
    }

    private void enter(@NotNull final GenomeRegion region) {
        start = region.start();
        end = region.end();
        ratioSupport = true;
    }

    private void exit(@NotNull final GenomeRegion region) {
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
            svStart = SegmentSupport.fromVariant(variant.type());
        }
    }

    private SegmentSupport svSupport(SegmentSupport current, StructuralVariantPosition variant) {
        return current.equals(NONE) ? SegmentSupport.fromVariant(variant.type()) : MULTIPLE;
    }

    private void insert(long end) {
        segments.add(ImmutablePurpleSegment.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .ratioSupport(ratioSupport)
                .support(svStart)
                .svCluster(false)
                .build());
    }
}
