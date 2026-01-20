package com.hartwig.hmftools.purple.segment;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

public class SegmentRefiner
{
    private final List<GenomeRegion> mSegments;

    public SegmentRefiner(final List<GenomeRegion> segments)
    {
        checkArgument(!segments.isEmpty());
        checkArgument(checkAllOnSameChromosome(segments), "Segments must be on the same chromosome");
        checkArgument(checkAreInOrder(segments), "Segments must be non-overlapping");
        this.mSegments = segments;
    }

    public List<GenomeRegion> segments()
    {
        return mSegments;
    }

    public List<PurpleSupportSegment> refine(List<PurpleSupportSegment> input)
    {
        checkArgument(!input.isEmpty());
        checkArgument(checkAllOnSameChromosome(input));
        checkArgument(Objects.equals(input.get(0).chromosome(), mSegments.get(0).chromosome()));
        checkArgument(checkAreInOrder(input));

        List<PurpleSupportSegment> result = new ArrayList<>();
        List<PurpleSupportSegment> remainingToSplit = new ArrayList<>(input);

        for(GenomeRegion segment : mSegments)
        {
            result.clear();
            for(PurpleSupportSegment pss : remainingToSplit)
            {
                result.addAll(pss.splitBy(segment));
            }
            remainingToSplit = new ArrayList<>(result);
        }
        return result;
    }

    private static <T extends GenomeRegion> boolean checkAllOnSameChromosome(final List<T> segments)
    {
        return segments.stream().map(GenomeRegion::chromosome).distinct().count() == 1;
    }

    private static <T extends GenomeRegion> boolean checkAreInOrder(List<T> segments)
    {
        for(int i = 1; i < segments.size(); ++i)
        {
            if(segments.get(i).start() <= segments.get(i - 1).end())
            {
                return false;
            }
        }
        return true;
    }
}
