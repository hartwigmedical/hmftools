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

        Map<String,List<ComparableItem>> sourceBreakends = Maps.newHashMap();

        for(Map.Entry<String,DatabaseAccess> entry : mConfig.DbConnections.entrySet())
        {
            sourceBreakends.put(entry.getKey(), getDisruptions(sampleId, entry.getValue()));
        }

        for(Map.Entry<String,List<ComparableItem>> entry1 : sourceBreakends.entrySet())
        {
            for(Map.Entry<String,List<ComparableItem>> entry2 : sourceBreakends.entrySet())
            {
                if(entry1 == entry2)
                    continue;

                final String source1 = entry1.getKey();
                final String source2 = entry2.getKey();

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, entry1.getValue(), entry2.getValue());
            }
        }
    }

    private List<ComparableItem> getDisruptions(final String sampleId, final DatabaseAccess dbAccess)
    {
        return dbAccess.readBreakends(sampleId).stream()
                .map(x -> new DisruptionData(x)).collect(Collectors.toList());
    }
}
