package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.List;
import java.util.Map;

public class SagaMatcher
{
    private final SagaResource mSagaResource;

    public SagaMatcher(final SagaResource sagaResource)
    {
        mSagaResource = sagaResource;
    }

    public SagaResource.Variant matchByLocation(final String chromosome, int position)
    {
        Map<String, List<SagaResource.IndexedBreakend>> breakends = mSagaResource.searchableBreakends();
        List<SagaResource.IndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchByLocationOnChromosome(position, chrBreakends);
    }

    public SagaResource.Variant matchByLocationOnChromosome(int position, List<SagaResource.IndexedBreakend> chrBreakends)
    {
        String bestVariant = null;
        int bestDistance = SAGA_LOCATION_MATCH_DISTANCE + 1;
        if(chrBreakends != null)
        {
            // TODO: can binary search for better performance
            for(SagaResource.IndexedBreakend breakend : chrBreakends)
            {
                int distance = abs(breakend.position().Position - position);
                if(distance < bestDistance)
                {
                    bestVariant = breakend.variantId();
                    bestDistance = distance;
                }
            }
        }
        if(bestVariant == null)
        {
            return null;

        }
        else
        {
            return mSagaResource.getVariantById(bestVariant);
        }
    }
}
