package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_RELATIVE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;

import java.util.List;

import com.google.common.collect.Lists;

public final class TincConstants
{
    protected static final int TINC_MIN_VARIANTS = 200;

    protected static final int TINC_GERMLINE_ABQ_MIN = 30;
    protected static final int TINC_MAX_FITTING_VARIANTS = 15_000;

    protected static final double TINC_GERMLINE_DEPTH_LOW = 0.5;
    protected static final double TINC_GERMLINE_DEPTH_HIGH = 1 + TINC_GERMLINE_DEPTH_LOW;

    protected static final int TINC_MQF_LIMIT = 25;

    protected static final int TINC_TUMOR_AF_LOWER_LIMIT = 50;
    protected static final int TINC_TUMOR_AF_UPPER_LIMIT = 300;
    protected static final double TINC_TUMOR_AF_UPPER_TEST_MIN = 0.011;

    public static final int TINC_GERMLINE_MAX_AD = 10;

    protected static final double TINC_RECOVERY_FACTOR = 2.5;
    protected static final double TINC_RECOVERY_MIN = 0.03;

    protected static final double TINC_RECOVERY_GERMLINE_AF_PROB = 0.001;

    protected static final List<String> RECOVERY_FILTERS = Lists.newArrayList(
            MAX_GERMLINE_RELATIVE_QUAL.filterName(), MAX_GERMLINE_VAF.filterName(), MAX_GERMLINE_ALT_SUPPORT.filterName(), PASS_FILTER);
}
