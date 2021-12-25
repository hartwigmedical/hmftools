package com.hartwig.hmftools.isofox.expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MULTI_MAP_QUALITY_THRESHOLD;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.calcGcRatioFromReadRegions;
import static com.hartwig.hmftools.isofox.common.FragmentType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;

public class ExpressionReadTracker
{
    private final IsofoxConfig mConfig;
    private final boolean mEnabled;
    private final GcRatioCounts mGcRatioCounts;
    private final GcRatioCounts mGeneGcRatioCounts;
    private final List<CategoryCountsData> mTransComboData;

    private GeneCollection mGenes;

    public ExpressionReadTracker(final IsofoxConfig config)
    {
        mConfig = config;
        mEnabled = config.runFunction(TRANSCRIPT_COUNTS);

        mGcRatioCounts = mConfig.requireGcRatioCalcs() ? new GcRatioCounts() : null;
        mGeneGcRatioCounts = mConfig.requireGcRatioCalcs() ? new GcRatioCounts() : null;
        mTransComboData = Lists.newArrayList();
        mGenes = null;
    }

    public boolean enabled() { return mEnabled; }
    public final GcRatioCounts getGcRatioCounts() { return mGcRatioCounts; }
    public final GcRatioCounts getGeneGcRatioCounts() { return mGeneGcRatioCounts; }
    public List<CategoryCountsData> getTransComboData() { return mTransComboData; }

    public void setGeneData(final GeneCollection genes)
    {
        if(mGeneGcRatioCounts != null)
            mGeneGcRatioCounts.clearCounts();

        mTransComboData.clear();

        mGenes = genes;
    }

    public void processUnsplicedGenes(
            final List<GeneReadData> overlapGenes, final List<Integer> validTranscripts, final List<int[]> commonMappings, int minMapQuality)
    {
        if(!mEnabled)
            return;

        List<String> unsplicedGeneIds = overlapGenes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        if(!unsplicedGeneIds.isEmpty())
        {
            CategoryCountsData catCounts = getCategoryCountsData(validTranscripts, unsplicedGeneIds);
            addGcCounts(catCounts, commonMappings, minMapQuality);
        }
    }

    public void processUnsplicedGenes(
            final FragmentMatchType comboTransMatchType, final List<GeneReadData> overlapGenes, final List<Integer> validTranscripts,
            final List<int[]> commonMappings, int minMapQuality)
    {
        if(!mEnabled)
            return;

        List<String> unsplicedGeneIds = comboTransMatchType == FragmentMatchType.SHORT ?
                overlapGenes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList()) : Lists.newArrayList();

        CategoryCountsData catCounts = getCategoryCountsData(validTranscripts, unsplicedGeneIds);
        addGcCounts(catCounts, commonMappings, minMapQuality);
    }

    public void processIntronicReads(final List<GeneReadData> genes, final ReadRecord read1, final ReadRecord read2)
    {
        if(!mEnabled)
            return;

        List<String> unsplicedGeneIds = genes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        if(!unsplicedGeneIds.isEmpty())
        {
            CategoryCountsData catCounts = getCategoryCountsData(Lists.newArrayList(), unsplicedGeneIds);

            List<int[]> readRegions = deriveCommonRegions(read1.getMappedRegionCoords(), read2.getMappedRegionCoords());
            addGcCounts(catCounts, readRegions, min(read1.mapQuality(), read2.mapQuality()));
        }
    }

    public void processEnrichedGeneFragments(int enrichedGeneFragments)
    {
        if(!mEnabled)
            return;

        // add to category counts
        final int[] enrichedRegion = mGenes.getEnrichedRegion();
        final List<String> unsplicedGeneIds = mGenes.findGenesCoveringRange(enrichedRegion[SE_START], enrichedRegion[SE_END], true)
                .stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        final List<Integer> transIds = mGenes.getEnrichedTranscripts().stream().map(x -> Integer.valueOf(x.TransId)).collect(Collectors.toList());
        CategoryCountsData catCounts = getCategoryCountsData(transIds, unsplicedGeneIds);

        // compute and cache GC data
        double gcRatio = calcGcRatioFromReadRegions(mConfig.RefGenome, mGenes.chromosome(), Lists.newArrayList(mGenes.getEnrichedRegion()));

        int[] gcRatioIndices = { -1, -1 };
        double[] gcRatioCounts = { 0, 0 };

        if(mGcRatioCounts != null)
        {
            mGcRatioCounts.determineRatioData(gcRatio, gcRatioIndices, gcRatioCounts);
            gcRatioCounts[0] *= enrichedGeneFragments;
            gcRatioCounts[1] *= enrichedGeneFragments;
        }

        addGcCounts(catCounts, gcRatioIndices, gcRatioCounts, enrichedGeneFragments);

    }

    private CategoryCountsData getCategoryCountsData(final List<Integer> transcripts, final List<String> geneIds)
    {
        CategoryCountsData transComboCounts = mTransComboData.stream()
                .filter(x -> x.matches(transcripts, geneIds)).findFirst().orElse(null);

        if(transComboCounts == null)
        {
            transComboCounts = new CategoryCountsData(transcripts, geneIds);

            if(mGcRatioCounts != null && mConfig.ApplyGcBiasAdjust)
                transComboCounts.initialiseGcRatioCounts(mGcRatioCounts.getCounts().length);

            mTransComboData.add(transComboCounts);
        }

        return transComboCounts;
    }


    public void processValidTranscript(int transId, final List<ReadRecord> reads, boolean isUniqueTrans)
    {
        final List<RegionReadData> processedRegions = Lists.newArrayList();

        for(ReadRecord read : reads)
        {
            List<RegionReadData> regions = read.getMappedRegions().entrySet().stream()
                    .filter(x -> x.getKey().hasTransId(transId))
                    .filter(x -> validExonMatch(x.getValue()))
                    .map(x -> x.getKey()).collect(Collectors.toList());

            for(RegionReadData region : regions)
            {
                if(!processedRegions.contains(region))
                {
                    // register a read against this valid transcript region
                    region.addTranscriptReadMatch(transId, isUniqueTrans);
                }
            }

            // any adjacent reads can record a splice junction count
            if(regions.size() > 1 && read.getTranscriptClassification(transId) == SPLICE_JUNCTION)
            {
                for(int r1 = 0; r1 < regions.size() - 1; ++r1)
                {
                    RegionReadData region1 = regions.get(r1);

                    for(int r2 = r1 + 1; r2 < regions.size(); ++r2)
                    {
                        RegionReadData region2 = regions.get(r2);

                        if(processedRegions.contains(region1) && processedRegions.contains(region2))
                            continue;

                        if(region1.getPostRegions().contains(region2))
                        {
                            region1.addTranscriptJunctionMatch(transId, SE_END, isUniqueTrans);
                            region2.addTranscriptJunctionMatch(transId, SE_START, isUniqueTrans);
                        }
                        else if(region1.getPreRegions().contains(region2))
                        {
                            region1.addTranscriptJunctionMatch(transId, SE_START, isUniqueTrans);
                            region2.addTranscriptJunctionMatch(transId, SE_END, isUniqueTrans);
                        }
                    }
                }
            }

            regions.forEach(x -> processedRegions.add(x));
        }
    }


    public void addGcCounts(final CategoryCountsData catCounts, final List<int[]> readRegions, int minMapQuality)
    {
        int[] gcRatioIndices = { -1, -1 };
        double[] gcRatioCounts = { 0, 0 };

        if (mGcRatioCounts != null)
        {
            double gcRatio = calcGcRatioFromReadRegions(mConfig.RefGenome, mGenes.chromosome(), readRegions);
            mGcRatioCounts.determineRatioData(gcRatio, gcRatioIndices, gcRatioCounts);
        }

        double fragmentCount = 1;

        if(minMapQuality <= MULTI_MAP_QUALITY_THRESHOLD && mConfig.ApplyMapQualityAdjust)
        {
            if(minMapQuality == 3)
                fragmentCount = 0.5;
            else if(minMapQuality == 2)
                fragmentCount = 0.33;
            else if(minMapQuality == 1)
                fragmentCount = 0.2;
            else
                fragmentCount = 0.1;

            mGenes.addCount(LOW_MAP_QUAL, 1);
        }

        addGcCounts(catCounts, gcRatioIndices, gcRatioCounts, fragmentCount);
    }

    public void addGcCounts(final CategoryCountsData catCounts, final int[] gcRatioIndices, double[] gcRatioCounts, double count)
    {
        if(mGcRatioCounts != null)
        {
            for(int i = 0; i < gcRatioIndices.length; ++i)
            {
                if (gcRatioIndices[i] >= 0)
                {
                    mGcRatioCounts.addGcRatioCount(gcRatioIndices[i], gcRatioCounts[i]);
                    mGeneGcRatioCounts.addGcRatioCount(gcRatioIndices[i], gcRatioCounts[i]);
                }
            }

            if(mConfig.ApplyGcBiasAdjust)
                catCounts.addGcRatioCounts(count, gcRatioIndices, gcRatioCounts);
            else
                catCounts.addCounts(count);
        }
        else
        {
            catCounts.addCounts(count);
        }
    }

}
