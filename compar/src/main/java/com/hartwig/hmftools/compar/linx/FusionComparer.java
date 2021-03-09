package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.FUSION;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
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

        Map<String,List<ComparableItem>> sourceFusions = Maps.newHashMap();

        for(Map.Entry<String,DatabaseAccess> entry : mConfig.DbConnections.entrySet())
        {
            sourceFusions.put(entry.getKey(), getFusions(sampleId, entry.getValue()));
        }

        for(Map.Entry<String,List<ComparableItem>> entry1 : sourceFusions.entrySet())
        {
            for(Map.Entry<String,List<ComparableItem>> entry2 : sourceFusions.entrySet())
            {
                if(entry1 == entry2)
                    continue;

                final String source1 = entry1.getKey();
                final String source2 = entry2.getKey();

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, entry1.getValue(), entry2.getValue());
            }
        }
    }

    private List<ComparableItem> getFusions(final String sampleId, final DatabaseAccess dbAccess)
    {
        return dbAccess.readFusions(sampleId).stream().map(x -> new FusionData(x)).collect(Collectors.toList());
    }


}
