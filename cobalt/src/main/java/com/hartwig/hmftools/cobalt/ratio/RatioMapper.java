package com.hartwig.hmftools.cobalt.ratio;

import tech.tablesaw.api.Table;

public interface RatioMapper
{
    // ratio mapper maps input ratio to output
    Table mapRatios(final Table inputRatios);
}
