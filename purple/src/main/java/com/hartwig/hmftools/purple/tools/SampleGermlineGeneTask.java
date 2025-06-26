package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.writeGeneDeletionData;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.writeGeneDeletionDetails;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.PurpleSegment;

public class SampleGermlineGeneTask implements Callable
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
    public Long call()
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

        return (long)0;
    }

    private void processSample(final String sampleId)
    {
        final Map<String,List<PurpleCopyNumber>> copyNumberMap = Maps.newHashMap();
        final Map<String,List<ObservedRegion>> fittedRegionMap = Maps.newHashMap();

        loadPurpleData(sampleId, copyNumberMap, fittedRegionMap);

        for(Map.Entry<String,List<ObservedRegion>> entry : fittedRegionMap.entrySet())
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
                PurpleCopyNumber matchedCopyNumber = copyNumbers.stream()
                        .filter(x -> positionsOverlap((int)x.start(), (int)x.end(), (int)region.start(), (int)region.end()))
                        .findFirst().orElse(null);

                if(matchedCopyNumber == null)
                {
                    PPL_LOGGER.error("sample({}) region({}: {}-{}) has no matching purple copy number",
                            sampleId, region.chromosome(), region.start(), region.end());
                    continue;
                }

                int regionLowerPos = region.start() - 500;
                int regionHighPos = region.end() + 500;

                // now find genes
                List<GeneData> overlappingGenes = geneDataList.stream()
                        .filter(x -> positionsOverlap(x.GeneStart, x.GeneEnd, regionLowerPos, regionHighPos))
                        .collect(Collectors.toList());

                StringJoiner geneNames = new StringJoiner(";");

                for(GeneData geneData : overlappingGenes)
                {
                    TranscriptData transData = mGeneDataCache.getTranscriptData(geneData.GeneId, "");

                    if(transData == null)
                        continue;

                    List<ExonData> overlappedExons = transData.exons().stream()
                            .filter(x -> positionsOverlap(x.Start, x.End, regionLowerPos, regionHighPos))
                            .collect(Collectors.toList());

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
            final String sampleId, final Map<String,List<PurpleCopyNumber>> copyNumberMap,
            final Map<String,List<ObservedRegion>> fittedRegionMap)
    {
        String samplePurpleDir = mPurpleDir.contains("*") ? mPurpleDir.replaceAll("\\*", sampleId) : mPurpleDir;

        try
        {
            List<PurpleCopyNumber> allCopyNumbers = PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(samplePurpleDir, sampleId));

            List<PurpleSegment> segments = PurpleSegment.read(PurpleSegment.generateFilename(samplePurpleDir, sampleId)).stream()
                .filter(x -> x.GermlineState == HET_DELETION || x.GermlineState == HOM_DELETION)
                .collect(Collectors.toList());

            for(PurpleSegment segment : segments)
            {
                List<ObservedRegion> regions = fittedRegionMap.get(segment.Chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    fittedRegionMap.put(segment.Chromosome, regions);

                    copyNumberMap.put(
                            segment.Chromosome,
                            allCopyNumbers.stream().filter(x -> x.chromosome().equals(segment.Chromosome)).collect(Collectors.toList()));
                }

                regions.add(ObservedRegion.fromSegment(segment));
            }

            PPL_LOGGER.debug("sample({}) read {} het-hom deletion regions",
                    sampleId, fittedRegionMap.values().stream().mapToInt(x -> x.size()).sum());
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("sample({}) failed to load purple files form {}: {}", sampleId, samplePurpleDir, e.toString());
        }
    }
}
