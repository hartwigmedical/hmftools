package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.LINX_DATA;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class LinxSvComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public LinxSvComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(LINX_DATA);

        Map<String,List<ComparableItem>> sourceSvData = Maps.newHashMap();

        for(Map.Entry<String, DatabaseAccess> entry : mConfig.DbConnections.entrySet())
        {
            sourceSvData.put(entry.getKey(), getSvData(sampleId, entry.getValue()));
        }

        for(Map.Entry<String,List<ComparableItem>> entry1 : sourceSvData.entrySet())
        {
            for(Map.Entry<String,List<ComparableItem>> entry2 : sourceSvData.entrySet())
            {
                if(entry1 == entry2)
                    continue;

                final String source1 = entry1.getKey();
                final String source2 = entry2.getKey();

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, entry1.getValue(), entry2.getValue());
            }
        }
    }

    private List<ComparableItem> getSvData(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<ComparableItem> svItems = Lists.newArrayList();

        final List<LinxCluster> clusters = dbAccess.readClusters(sampleId);
        final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
        final List<LinxSvAnnotation> annotations = dbAccess.readSvAnnotations(sampleId);

        for(StructuralVariantData svData : svDataList)
        {
            LinxSvAnnotation annotation = annotations.stream().filter(x -> x.svId() == svData.id()).findFirst().orElse(null);

            if(annotation == null)
                continue;

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

            if(cluster == null)
                continue;

            svItems.add(new LinxSvData(svData, annotation, cluster));
        }

        return svItems;
    }


}
