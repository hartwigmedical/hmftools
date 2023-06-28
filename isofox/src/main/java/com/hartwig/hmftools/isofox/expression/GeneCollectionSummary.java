package com.hartwig.hmftools.isofox.expression;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.results.GeneResult;
import com.hartwig.hmftools.isofox.results.TranscriptResult;

public class GeneCollectionSummary
{
    public final String ChrId;
    public final List<String> GeneIds;
    public final String GeneNames;
    public final List<CategoryCountsData> TransCategoryCounts;

    public final List<GeneResult> GeneResults;
    public final List<TranscriptResult> TranscriptResults;

    private final int mTotalLowMqFragments;
    private final Map<String,Double> mFitAllocations; // results from the expected rate vs counts fit routine, stored per transcript
    private double mFitResiduals;

    public GeneCollectionSummary(
            final String chrId, final List<String> geneIds, final String geneNames, final List<CategoryCountsData> transCategoryCounts)
    {
        ChrId = chrId;
        GeneIds = geneIds;
        GeneNames = geneNames;
        TransCategoryCounts = Lists.newArrayList(transCategoryCounts);
        GeneResults = Lists.newArrayList();
        TranscriptResults = Lists.newArrayList();

        mTotalLowMqFragments = TransCategoryCounts.stream().mapToInt(x -> x.lowMapQualFragments()).sum();

        mFitAllocations = Maps.newHashMap();
        mFitResiduals = 0;
    }

    public int spliceGenesCount() { return (int)GeneResults.stream().filter(x -> x.getSplicedAlloc() > 0).count(); }

    public void setFitResiduals(double residuals) { mFitResiduals = residuals; }
    public double getFitResiduals() { return mFitResiduals; }

    public Map<String,Double> getFitAllocations() { return mFitAllocations; }

    public double getFitAllocation(final String transName)
    {
        Double allocation = mFitAllocations.get(transName);
        return allocation != null ? allocation : 0;
    }

    public void setFitAllocations()
    {
        Map<String,Double> geneSpliceTotals = Maps.newHashMap();

        for (final TranscriptResult transResult : TranscriptResults)
        {
            final String transName = transResult.Trans.TransName;
            double fitAllocation = getFitAllocation(transName);

            transResult.setFitAllocation(fitAllocation);

            Double geneFitAllocation = geneSpliceTotals.get(transResult.Trans.GeneId);
            if(geneFitAllocation == null)
                geneSpliceTotals.put(transResult.Trans.GeneId, fitAllocation);
            else
                geneSpliceTotals.put(transResult.Trans.GeneId, geneFitAllocation + fitAllocation);

        }

        for(final GeneResult geneResult : GeneResults)
        {
            Double geneFitAllocation = geneSpliceTotals.get(geneResult.Gene.GeneId);
            geneResult.setFitAllocation(
                    geneFitAllocation != null ? geneFitAllocation : 0, getFitAllocation(geneResult.Gene.GeneId));
        }
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
            long totalGeneLength = GeneResults.stream().mapToLong(x -> x.Gene.length()).sum();

            for (final GeneResult geneResult : GeneResults)
            {
                double residualsFraction = geneResult.Gene.length() / (double) totalGeneLength * mFitResiduals;
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
        // eg if old total was 10K and new is 2K, then will multiply all new counts by 5
        double adjustFactor = originalTotal/newTotal;
        TransCategoryCounts.forEach(x -> x.adjustCounts(adjustFactor));
    }

    public void assignLowMapQualityFragments()
    {
        if(mTotalLowMqFragments == 0)
            return;

        double totalTranscriptAlloc = TranscriptResults.stream().mapToDouble(x -> x.getFitAllocation()).sum();
        double totalUnsplicedAlloc = GeneResults.stream().mapToDouble(x -> x.getUnsplicedAlloc()).sum();
        double totalAlloc = totalTranscriptAlloc + totalUnsplicedAlloc;

        if(totalAlloc == 0)
            return;

        double splicedLowMqFrags = mTotalLowMqFragments * totalTranscriptAlloc / totalAlloc;
        double unsplicedLowMqFrags = mTotalLowMqFragments * totalUnsplicedAlloc / totalAlloc;

        // divide amongst transcripts
        for(final TranscriptResult transResult : TranscriptResults)
        {
            double transAlloc = totalTranscriptAlloc > 0 ? transResult.getFitAllocation() / totalTranscriptAlloc * splicedLowMqFrags : 0;
            transResult.setLowMapQualsAllocation(transAlloc);
        }

        // split amongst genes as per fragment allocation
        for(final GeneResult geneResult : GeneResults)
        {
            double splicedAlloc = totalTranscriptAlloc > 0 ? geneResult.getSplicedAlloc() / totalTranscriptAlloc * splicedLowMqFrags : 0;
            double unsplicedAlloc = totalUnsplicedAlloc > 0 ? geneResult.getUnsplicedAlloc() / totalUnsplicedAlloc * unsplicedLowMqFrags : 0;
            geneResult.setLowMapQualsAllocation(splicedAlloc + unsplicedAlloc);
        }
    }
}
