package com.hartwig.hmftools.purple;

import static java.lang.Math.max;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.purple.segment.PurpleSegmentFactory.validateSegments;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.region.ObservedRegionFactory;
import com.hartwig.hmftools.purple.segment.PurpleSegment;
import com.hartwig.hmftools.purple.segment.PurpleSegmentFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.purple.segment.PCFPositionsSupplier;

public class Segmentation
{
    private final Multimap<Chromosome, GCProfile> mGcProfiles;
    private final ReferenceData mReferenceData;
    private final int mWindowSize;

    public Segmentation(final ReferenceData referenceData) throws IOException
    {
        mReferenceData = referenceData;
        mWindowSize = WINDOW_SIZE;

        PPL_LOGGER.info("reading GC Profiles from {}", referenceData.GcProfileFilename);
        mGcProfiles = GCProfileFactory.loadGCContent(mWindowSize, referenceData.GcProfileFilename);
    }

    public List<ObservedRegion> createObservedRegions(
            final List<StructuralVariant> structuralVariants, final AmberData amberData, final CobaltData cobaltData)
    {
        final Map<Chromosome,List<PCFPosition>> pcfPositions = PCFPositionsSupplier.createPositions(amberData, cobaltData);

        final PurpleSegmentFactory factory = new PurpleSegmentFactory(
                mWindowSize, mReferenceData.Centromeres, mReferenceData.ChromosomeLengths);

        final List<PurpleSegment> segments = factory.createSegments(structuralVariants, pcfPositions, cobaltData.Ratios);

        if(!validateSegments(segments))
            return Lists.newArrayList();

        final ObservedRegionFactory observedRegionFactory = new ObservedRegionFactory(mWindowSize, cobaltData.CobaltChromosomes);

        List<ObservedRegion> observedRegions = observedRegionFactory.formObservedRegions(
                segments, amberData.ChromosomeBafs, cobaltData.Ratios, mGcProfiles);

        return observedRegions;
    }

    public static boolean validateObservedRegions(final List<ObservedRegion> observedRegions)
    {
        boolean isValid = true;

        for(int i = 1; i < observedRegions.size(); ++i)
        {
            ObservedRegion observedRegion = observedRegions.get(i);

            if(observedRegion.start() > observedRegion.end())
            {
                PPL_LOGGER.error("observed region({}:{}-{}) is invalid",
                        observedRegion.chromosome(), observedRegion.start(), observedRegion.end());

                isValid = false;
            }

            /* the start pos can fall within (at the mid-point to be precise) the min and max start
            else if(!positionsWithin(observedRegion.minStart(), observedRegion.maxStart(), observedRegion.start(), observedRegion.end()))
            {
                PPL_LOGGER.error("observed region({}:{}-{}) has invalid min/maxStart({}-{})",
                        observedRegion.chromosome(), observedRegion.start(), observedRegion.end(),
                        observedRegion.minStart(), observedRegion.maxStart());

                isValid = false;
            }
            */

            ObservedRegion prevRegion = observedRegions.get(i - 1);

            if(!observedRegion.chromosome().equals(prevRegion.chromosome()))
                continue;

            if(observedRegion.start() <= prevRegion.end())
            {
                PPL_LOGGER.error("observed region({}:{}-{}) overlaps previous({}-{})",
                        observedRegion.chromosome(), observedRegion.start(), observedRegion.end(), prevRegion.start(), prevRegion.end());
                isValid = false;
            }
        }

        return isValid;
    }
}
