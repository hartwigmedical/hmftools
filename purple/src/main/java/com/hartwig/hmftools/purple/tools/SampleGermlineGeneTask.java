package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findMatchingCopyNumber;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findOverlappingExons;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findOverlappingGenes;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.loadPurpleDataForGermline;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.writeGeneDeletionData;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.writeGeneDeletionDetails;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class SampleGermlineGeneTask implements Callable<Void>
{
    private final int mTaskId;
    private final BufferedWriter mWriter;
    private final EnsemblDataCache mGeneDataCache;
    private final String mPurpleDir;

    private final List<String> mSampleIds;
    private final boolean mWriteVerbose;

    public SampleGermlineGeneTask(
            int taskId, final BufferedWriter writer, final EnsemblDataCache geneDataCache, final String purpleDir, final boolean writeVerbose)
    {
        mTaskId = taskId;
        mGeneDataCache = geneDataCache;
        mWriter = writer;
        mPurpleDir = purpleDir;
        mWriteVerbose = writeVerbose;

        mSampleIds = Lists.newArrayList();
    }

    public List<String> getSampleIds() { return mSampleIds; }

    @Override
    public Void call()
    {
        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            processSample(sampleId);

            if(i > 0 && (i % 100) == 0)
            {
                PPL_LOGGER.debug("{}: processed {} samples", mTaskId, i);
            }
        }

        PPL_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());

        return null;
    }

    private void processSample(final String sampleId)
    {
        final Map<String, List<PurpleCopyNumber>> copyNumberMap = Maps.newHashMap();
        final Map<String, List<ObservedRegion>> fittedRegionMap = Maps.newHashMap();

        loadPurpleData(sampleId, copyNumberMap, fittedRegionMap);

        for(Map.Entry<String, List<ObservedRegion>> entry : fittedRegionMap.entrySet())
        {
            String chromosome = entry.getKey();

            List<PurpleCopyNumber> copyNumbers = copyNumberMap.get(chromosome);
            List<GeneData> geneDataList = mGeneDataCache.getChrGeneDataMap().get(chromosome);

            if(copyNumbers == null || geneDataList == null)
                continue;

            for(ObservedRegion region : entry.getValue())
            {
                if(region.germlineStatus() != HOM_DELETION && region.germlineStatus() != HET_DELETION)
                    continue;

                // find the overlapping / matching copy number region
                PurpleCopyNumber matchedCopyNumber = findMatchingCopyNumber(region, copyNumbers);

                if(matchedCopyNumber == null)
                {
                    PPL_LOGGER.error("sample({}) region({}: {}-{}) has no matching purple copy number",
                            sampleId, region.chromosome(), region.start(), region.end());
                    continue;
                }

                final int BUFFER = 500;

                // now find genes
                List<GeneData> overlappingGenes = findOverlappingGenes(
                        region.chromosome(), region.start(), region.end(), BUFFER, geneDataList);

                StringJoiner geneNames = new StringJoiner(";");

                for(GeneData geneData : overlappingGenes)
                {
                    TranscriptData transData = mGeneDataCache.getCanonicalTranscriptData(geneData.GeneId);

                    if(transData == null)
                        continue;

                    List<ExonData> overlappedExons = findOverlappingExons(
                            transData, region.start(), region.end(), BUFFER);

                    if(overlappedExons.isEmpty())
                        continue;

                    geneNames.add(geneData.GeneName);

                    PPL_LOGGER.trace("sample({}) region({}: {}-{}) overlaps gene({}) exons({})",
                            sampleId, region.chromosome(), region.start(), region.end(), geneData.GeneName, overlappedExons.size());

                    if(mWriteVerbose)
                    {
                        writeGeneDeletionDetails(mWriter, sampleId, region, matchedCopyNumber, geneData, transData, overlappedExons);
                    }
                }

                if(!mWriteVerbose)
                {
                    writeGeneDeletionData(mWriter, sampleId, region, geneNames.toString());
                }
            }
        }
    }

    private void loadPurpleData(
            final String sampleId, final Map<String, List<PurpleCopyNumber>> copyNumberMap,
            final Map<String, List<ObservedRegion>> fittedRegionMap)
    {
        // Call the utility method with all required parameters
        loadPurpleDataForGermline(sampleId, mPurpleDir, copyNumberMap, fittedRegionMap);
    }
}

class GeneAmplification
{
    private final int totalExonCount;
    private final int overlappingExonCount;
    private final boolean overlapsFirstExon;
    private final boolean overlapsLastExon;

    GeneAmplification(List<ExonData> allExons, int start, int end)
    {
        Preconditions.checkArgument(!allExons.isEmpty());
        Preconditions.checkArgument(start < end);

        totalExonCount = allExons.size();
        List<ExonData> overlappingExons = allExons.stream().filter(x -> positionsOverlap(x.Start, x.End, start, end)).toList();
        overlappingExonCount = overlappingExons.size();

        if (overlappingExonCount > 0)
        {
            overlapsFirstExon = overlappingExons.get(0).equals(allExons.get(0));
            overlapsLastExon = overlappingExons.get(overlappingExonCount - 1).equals(allExons.get(totalExonCount - 1));
        }
        else
        {
            overlapsFirstExon = false;
            overlapsLastExon = false;
        }
    }

    public boolean isCompleteAmplification()
    {
        return overlappingExonCount == totalExonCount;
    }

    public boolean isOfInterest()
    {
        if (overlappingExonCount == 0)
        {
            return false;
        }
        if (isCompleteAmplification())
        {
            return true;
        }
        return !overlapsFirstExon && !overlapsLastExon;
    }

    public boolean isTailAmplification()
    {
        return overlappingExonCount > 0 && !isCompleteAmplification() && overlapsLastExon;
    }

    public boolean isHeadAmplification()
    {
        return overlappingExonCount > 0 && !isCompleteAmplification() && overlapsFirstExon;
    }

    public int numberOfAffectedExons()
    {
        return overlappingExonCount;
    }
}
