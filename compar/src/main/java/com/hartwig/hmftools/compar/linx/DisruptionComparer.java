package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.Category.DRIVER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class DisruptionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DisruptionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(DISRUPTION);

        final List<List<ComparableItem>> sourceBreakends = Lists.newArrayList();

        for(String sourceName : mConfig.DbSourceNames)
        {
            sourceBreakends.add(getDisruptions(sampleId, mConfig.DbConnections.get(sourceName)));
        }

        for(int i = 0; i < mConfig.DbSourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.DbSourceNames.get(i);

            for(int j = i + 1; j < mConfig.DbSourceNames.size(); ++j)
            {
                final String source2 = mConfig.DbSourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourceBreakends.get(i), sourceBreakends.get(j));
            }
        }
    }

    private List<ComparableItem> getDisruptions(final String sampleId, final DatabaseAccess dbAccess)
    {
        return dbAccess.readBreakends(sampleId).stream()
                .map(x -> new DisruptionData(x)).collect(Collectors.toList());
    }
}
