package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class FusionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public FusionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(FUSION);

        final List<List<ComparableItem>> sourceFusions = Lists.newArrayList();

        for(String sourceName : mConfig.DbSourceNames)
        {
            sourceFusions.add(getFusions(sampleId, mConfig.DbConnections.get(sourceName)));
        }

        for(int i = 0; i < mConfig.DbSourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.DbSourceNames.get(i);

            for(int j = i + 1; j < mConfig.DbSourceNames.size(); ++j)
            {
                final String source2 = mConfig.DbSourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourceFusions.get(i), sourceFusions.get(j));
            }
        }
    }

    private List<ComparableItem> getFusions(final String sampleId, final DatabaseAccess dbAccess)
    {
        return dbAccess.readFusions(sampleId).stream().map(x -> new FusionData(x)).collect(Collectors.toList());
    }


}
