package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.common.sigs.DataUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.gc.GcRatioCounts;
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

    private final Map<String,Double> mFitAllocations; // results from the expected rate vs counts fit routine, stored per transcript
    private final Map<String,Double> mMedianExpectedGcRatios; //
    private double mFitResiduals;

    public GeneCollectionSummaryData(
            final String chrId, final List<String> geneIds, final String geneNames, final List<CategoryCountsData> transCategoryCounts)
    {
        ChrId = chrId;
        GeneIds = geneIds;
        GeneNames = geneNames;
        TransCategoryCounts = Lists.newArrayList(transCategoryCounts);
        GeneResults = Lists.newArrayList();
        TranscriptResults = Lists.newArrayList();

        mFitAllocations = Maps.newHashMap();
        mMedianExpectedGcRatios = Maps.newHashMap();
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

    public void addMedianExpectedGcRatio(final String transName, double gcRatio)
    {
        mMedianExpectedGcRatios.put(transName, gcRatio);
    }

    public double getMedianExpectedGcRatio(final String transName)
    {
        Double gcRatio = mMedianExpectedGcRatios.get(transName);
        return gcRatio != null ? gcRatio : 0;
    }

    public void allocateResidualsToGenes()
    {
        if(GeneResults.size() == 1)
        {
            GeneResults.get(0).setFitResiduals(mFitResiduals);
        }
        else
        {
            // divvy up residuals between the genes according to their length
            long totalGeneLength = GeneResults.stream().mapToLong(x -> x.geneData().length()).sum();

            for (final GeneResult geneResult : GeneResults)
            {
                double residualsFraction = geneResult.geneData().length() / (double) totalGeneLength * mFitResiduals;
                geneResult.setFitResiduals(residualsFraction);
            }
        }
    }

    public void applyGcAdjustments(final double[] gcAdjustments)
    {
        double originalTotal = 0;
        double newTotal = 0;

        for(final CategoryCountsData catCounts : TransCategoryCounts)
        {
            originalTotal += catCounts.fragmentCount();
            catCounts.applyGcAdjustments(gcAdjustments);
            newTotal += catCounts.fragmentCount();
        }

        // ensure no overall net change to counts after the adjustment
        // eg if old total was 10K and new is 2K, then will mutiply all new counts by 5
        double adjustFactor = originalTotal/newTotal;
        TransCategoryCounts.forEach(x -> x.adjustCounts(adjustFactor));
    }

    public void recordMedianGcRatios()
    {
        for (final GeneResult geneResult : GeneResults)
        {
            GcRatioCounts gcCounts = new GcRatioCounts();

            for(final CategoryCountsData catCounts : TransCategoryCounts)
            {
                if(catCounts.unsplicedGeneIds().contains(geneResult.geneData().GeneId))
                    gcCounts.mergeRatioCounts(catCounts.fragmentCountsByGcRatio());
            }

            geneResult.setMedianActualGcRatio(gcCounts.getPercentileRatio(0.5));
        }
    }

}
