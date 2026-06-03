package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BasePosition;

import org.jetbrains.annotations.Nullable;

// Matches junctions to SAGA variants by genomic location.
// For a junction, find the nearest SAGA variant breakend within the configured distance tolerance.
public class SagaLocationMatcher
{
    private final Config mConfig;
    // Map from contig to list of breakends sorted by position ascending.
    private final Map<String, List<SagaIndexedBreakend>> mSearchableBreakends;

    // Do not instantiate directly - use SagaMatcherFactory.
    // searchableBreakends MUST be sorted ascending by breakend position.
    public SagaLocationMatcher(final Config config, final Map<String, List<SagaIndexedBreakend>> searchableBreakends)
    {
        mConfig = config;
        mSearchableBreakends = searchableBreakends;
    }

    public record Config(
            int locationDistanceMax
    )
    {
        public static final Config DEFAULT = new Config(SAGA_LOCATION_MATCH_DISTANCE);
    }

    @Nullable
    public SagaMatchByLocation match(final String chromosome, int position)
    {
        Map<String, List<SagaIndexedBreakend>> breakends = mSearchableBreakends;
        List<SagaIndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchOnChromosome(position, chrBreakends);
    }

    @Nullable
    private SagaMatchByLocation matchOnChromosome(int position, final List<SagaIndexedBreakend> chrBreakends)
    {
        SagaIndexedBreakend bestBreakend = null;
        int bestDistance = mConfig.locationDistanceMax + 1;
        if(chrBreakends != null && !chrBreakends.isEmpty())
        {
            SagaIndexedBreakend target =
                    new SagaIndexedBreakend(new BasePosition(chrBreakends.get(0).chromosome(), position), null);
            int index = Collections.binarySearch(chrBreakends, target);
            if(index >= 0)
            {
                // Exact position found.
                bestBreakend = chrBreakends.get(index);
                bestDistance = 0;
            }
            else
            {
                // Exact position not found, so need to look at adjacent indices.
                index = -(index + 1);
                int startIndex = max(index - 1, 0);
                int endIndex = min(index + 1, chrBreakends.size() - 1);
                for(int i = startIndex; i <= endIndex; i++)
                {
                    SagaIndexedBreakend breakend = chrBreakends.get(i);
                    int distance = abs(breakend.position() - position);
                    if(distance < bestDistance)
                    {
                        bestBreakend = breakend;
                        bestDistance = distance;
                    }
                }
            }
        }
        if(bestBreakend == null)
        {
            return null;
        }
        else
        {
            SagaVariant variant = bestBreakend.variant();
            return new SagaMatchByLocation(variant, bestDistance);
        }
    }
}
