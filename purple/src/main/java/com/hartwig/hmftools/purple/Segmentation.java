package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.region.ObservedRegionFactory;
import com.hartwig.hmftools.purple.segment.PurpleSegment;
import com.hartwig.hmftools.purple.segment.PurpleSegmentFactory;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.segment.PCFPositionsSupplier;

import org.jetbrains.annotations.NotNull;

class Segmentation
{
    private final Multimap<Chromosome, GCProfile> mGcProfiles;
    private final ReferenceData mReferenceData;
    private final int mWindowSize;

    public Segmentation(@NotNull final PurpleConfig config, final ReferenceData referenceData) throws IOException
    {
        mReferenceData = referenceData;
        mWindowSize = config.commonConfig().windowSize();

        PPL_LOGGER.info("Reading GC Profiles from {}", config.commonConfig().gcProfile());
        mGcProfiles = GCProfileFactory.loadGCContent(mWindowSize, config.commonConfig().gcProfile());
    }

    @NotNull
    public List<ObservedRegion> createSegments(
            final List<StructuralVariant> structuralVariants, final AmberData amberData, final CobaltData cobaltData)
    {
        final Multimap<Chromosome, PCFPosition> pcfPositions = PCFPositionsSupplier.createPositions(amberData, cobaltData);

        final PurpleSegmentFactory factory = new PurpleSegmentFactory(mWindowSize,
                mReferenceData.Centromeres, mReferenceData.ChromosomeLengths);

        final List<PurpleSegment> segments = factory.segment(structuralVariants, pcfPositions, cobaltData.Ratios);

        final ObservedRegionFactory observedRegionFactory =
                new ObservedRegionFactory(mWindowSize, cobaltData.CobaltChromosomes);

        return observedRegionFactory.combine(segments, amberData.ChromosomeBafs, cobaltData.Ratios, mGcProfiles);
    }
}
