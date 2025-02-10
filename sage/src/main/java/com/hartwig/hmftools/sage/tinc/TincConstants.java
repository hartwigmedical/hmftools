package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_RELATIVE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;

import java.util.List;

import com.google.common.collect.Lists;

public final class TincConstants
{
    protected static final int TINC_GERMLINE_ABQ_MIN = 30;
    protected static final int TINC_MAX_FITTING_VARIANTS = 15_000;

    protected static final double TINC_GERMLINE_DEPTH_LOW = 0.5;
    protected static final double TINC_GERMLINE_DEPTH_HIGH = 1 + TINC_GERMLINE_DEPTH_LOW;

    protected static final int TINC_GERMLINE_MAX_AD = 10;

    protected static final double TINC_RECOVERY_FACTOR = 2.5;
    protected static final double TINC_RECOVERY_MIN = 0.03;

    protected static final List<String> RECOVERY_FILTERS = Lists.newArrayList(
            MAX_GERMLINE_RELATIVE_QUAL.filterName(), MAX_GERMLINE_VAF.filterName(), MAX_GERMLINE_ALT_SUPPORT.filterName(), PASS);
}
