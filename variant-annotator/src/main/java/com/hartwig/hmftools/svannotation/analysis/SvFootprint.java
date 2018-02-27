package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svannotation.analysis.SvClusterData;

import java.util.List;

public class SvFootprint {

    private int mFootprintId;

    private List<SvClusterData> mClusteredSVs;

    public SvFootprint(final int footprintId)
    {
        mFootprintId = footprintId;
        mClusteredSVs = Lists.newArrayList();
    }

    public int getFootprintId() { return mFootprintId; }

    public List<SvClusterData> getSVs() { return mClusteredSVs; }

    public void addVariant(final SvClusterData variant)
    {
        mClusteredSVs.add(variant);
    }
}
