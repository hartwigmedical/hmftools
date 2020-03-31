package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.TranscriptResult;

public class GeneCollectionSummaryData
{
    public final String ChrId;
    public final List<String> GeneIds;
    public final String GeneNames;
    public final List<CategoryCountsData> TransCategoryCounts;

    public final List<GeneResult> GeneResults;
    public final List<TranscriptResult> TranscriptResults;

    private final Map<String,Double> mFitAllocations; // results from the expected rate vs counts fit routine
    private double mFitResiduals;

    public GeneCollectionSummaryData(
            final String chrId, final List<String> geneIds, final String geneNames, final List<CategoryCountsData> transCategoryCounts)
    {
        ChrId = chrId;
        GeneIds = geneIds;
        GeneNames = geneNames;
        TransCategoryCounts = transCategoryCounts;
        GeneResults = Lists.newArrayList();
        TranscriptResults = Lists.newArrayList();

        mFitAllocations = Maps.newHashMap();
        mFitResiduals = 0;
    }

    public void setFitResiduals(double residuals) { mFitResiduals = residuals; }
    public double getFitResiduals() { return mFitResiduals; }

    public Map<String,Double> getFitAllocations() { return mFitAllocations; }

    public double getFitAllocation(final String transName)
    {
        Double allocation = mFitAllocations.get(transName);
        return allocation != null ? allocation : 0;
    }

}
