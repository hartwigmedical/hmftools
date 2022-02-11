package com.hartwig.hmftools.sage.filters;

import static com.hartwig.hmftools.sage.config.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.config.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.config.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.config.SoftFilter.isGermlineAndNotTumorFiltered;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;

import org.junit.Test;

public class SoftFilterTest
{
    @Test
    public void testIsGermlineAndNotTumorFiltered()
    {
        assertTrue(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.toString())));
        assertTrue(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.toString(), MAX_GERMLINE_VAF.toString())));

        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet()));
        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.toString(), MIN_TUMOR_QUAL.toString())));
        assertFalse(isGermlineAndNotTumorFiltered(Sets.newHashSet(MAX_GERMLINE_REL_RAW_BASE_QUAL.toString(), "Random Filter")));
    }
}
