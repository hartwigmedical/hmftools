package com.hartwig.hmftools.sage.common;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.config.VariantFilters;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class SageVariantFactory
{
    private final VariantFilters mFilters;

    public SageVariantFactory(final FilterConfig config)
    {
        mFilters = new VariantFilters(config);
    }

    public SageVariant create(final Candidate candidate, final List<ReadContextCounter> normal, final List<ReadContextCounter> tumor)
    {
        boolean isNormalEmpty = normal.isEmpty();

        final VariantTier tier = candidate.tier();
        final SoftFilterConfig softFilterConfig = mFilters.getTieredSoftFilterConfig(tier);

        SageVariant variant = new SageVariant(candidate, normal, tumor);

        if(!mFilters.enabled())
        {
            return new SageVariant(candidate, normal, tumor);
        }

        final Set<String> variantFilters = variant.filters();

        // where there are multiple tumor samples, if any of them pass then clear any filters from the others
        for(ReadContextCounter tumorReadContextCounter : tumor)
        {
            final Set<String> tumorFilters = Sets.newHashSet();

            mFilters.applyTumorFilters(tier, softFilterConfig, tumorReadContextCounter, tumorFilters);

            if(!isNormalEmpty)
            {
                mFilters.applyTumorNormalFilters(tier, softFilterConfig, normal.get(0), tumorReadContextCounter, tumorFilters);
            }

            if(tumorFilters.isEmpty())
            {
                variantFilters.clear();
                break;
            }
            else
            {
                variantFilters.addAll(tumorFilters);
            }
        }

        return new SageVariant(candidate, normal, tumor);
    }
}
