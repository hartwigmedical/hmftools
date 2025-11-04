package com.hartwig.hmftools.qsee.prep;

import java.util.List;

import com.hartwig.hmftools.qsee.prep.category.BamMetricsPrep;
import com.hartwig.hmftools.qsee.prep.category.CobaltGcMediansPrep;
import com.hartwig.hmftools.qsee.prep.category.ReduxBqrPrep;
import com.hartwig.hmftools.qsee.prep.category.ReduxMsIndelErrorPrep;
import com.hartwig.hmftools.qsee.prep.category.SummaryTablePrep;

public class CategoryPrepFactory
{
    PrepConfig mConfig;

    public CategoryPrepFactory(PrepConfig config)
    {
        mConfig = config;
    }

    public List<CategoryPrep> createCategoryPreps()
    {
        return List.of(
                new SummaryTablePrep(mConfig),
                new BamMetricsPrep(mConfig),
                new CobaltGcMediansPrep(mConfig),
                new ReduxBqrPrep(mConfig),
                new ReduxMsIndelErrorPrep(mConfig)
        );
    }
}
