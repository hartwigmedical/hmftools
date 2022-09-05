package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;

import java.io.IOException;
import java.util.List;
import java.util.Map;

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

    public List<ObservedRegion> createSegments(
            final List<StructuralVariant> structuralVariants, final AmberData amberData, final CobaltData cobaltData)
    {
        final Map<Chromosome,List<PCFPosition>> pcfPositions = PCFPositionsSupplier.createPositions(amberData, cobaltData);

        final PurpleSegmentFactory factory = new PurpleSegmentFactory(
                mWindowSize, mReferenceData.Centromeres, mReferenceData.ChromosomeLengths);

        final List<PurpleSegment> segments = factory.segment(structuralVariants, pcfPositions, cobaltData.Ratios);

        final ObservedRegionFactory observedRegionFactory = new ObservedRegionFactory(mWindowSize, cobaltData.CobaltChromosomes);

        List<ObservedRegion> observedRegions = observedRegionFactory.combine(segments, amberData.ChromosomeBafs, cobaltData.Ratios, mGcProfiles);

        // correct the start position to the mid-point of the min and max start
        for(int i = 0; i < observedRegions.size(); ++i)
        {
            ObservedRegion observedRegion = observedRegions.get(i);

            if(observedRegion.minStart() < observedRegion.maxStart())
            {
                int midpoint = (observedRegion.minStart() + observedRegion.maxStart()) / 2;
                observedRegion.setStart(midpoint);

                if(i > 0)
                {
                    ObservedRegion prevRegion = observedRegions.get(i - 1);

                    if(prevRegion.chromosome().equals(observedRegion.chromosome()))
                    {
                        prevRegion.setEnd(midpoint - 1);
                    }
                }
            }
        }

        return observedRegions;
    }
}
