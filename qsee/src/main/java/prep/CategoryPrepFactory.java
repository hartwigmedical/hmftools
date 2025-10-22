package prep;

import java.util.List;

import prep.category.CobaltGcMediansPrep;
import prep.category.ReduxBqrPrep;

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
                new CobaltGcMediansPrep(mConfig),
                new ReduxBqrPrep(mConfig)
        );
    }
}
