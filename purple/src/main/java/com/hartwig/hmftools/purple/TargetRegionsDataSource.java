package com.hartwig.hmftools.purple;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsDataSource implements TargetRegionsCopyNumbers.DataSupplier
{
    private final TargetRegionsData mReferenceDataTargetRegions;
    private final RefGenomeVersion mRefGenomeVersion;
    private final List<PurpleSegment> mSegments;
    private PurpleSegment mLastSegmentReturned = null;

    public TargetRegionsDataSource(final TargetRegionsData targetRegionsData, final RefGenomeVersion genomeVersion,
            final List<PurpleSegment> mSegments)
    {
        this.mReferenceDataTargetRegions = targetRegionsData;
        this.mRefGenomeVersion = genomeVersion;
        this.mSegments = mSegments;
    }

    @Override
    public List<TaggedRegion> targetRegions(final Chromosome chromosome)
    {
        return mReferenceDataTargetRegions.targetRegions(mRefGenomeVersion.versionedChromosome(chromosome));
    }

    @Override
    public GermlineStatus germlineStatus(final GenomePosition position)
    {
        // Successive calls to this method are generally for successive positions
        // (and from a single thread) so caching is possible.
        if(mLastSegmentReturned != null && mLastSegmentReturned.containsPosition(position))
        {
            return mLastSegmentReturned.GermlineState;
        }
        for(PurpleSegment segment : mSegments)
        {
            if(segment.containsPosition(position))
            {
                mLastSegmentReturned = segment;
                return segment.GermlineState;
            }
        }
        return null;
    }
}
