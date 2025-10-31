package prep;

import java.util.List;

import prep.category.BamMetricsPrep;
import prep.category.CobaltGcMediansPrep;
import prep.category.ReduxBqrPrep;
import prep.category.ReduxMsIndelErrorPrep;
import prep.category.SummaryTablePrep;

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
