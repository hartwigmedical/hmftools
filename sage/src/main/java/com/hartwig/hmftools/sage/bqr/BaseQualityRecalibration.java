package com.hartwig.hmftools.sage.bqr;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.sage.SageCommon.generateBqrFilename;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BaseQualityRecalibration
{
    private final SageConfig mConfig;
    private final String mPanelBedFile;
    private final List<String> mTumorIds;
    private final List<String> mTumorBams;
    private final IndexedFastaSequenceFile mRefGenome;

    private final Map<String, BqrRecordMap> mSampleRecalibrationMap;
    private final Queue<PartitionTask> mRegions;
    private final BaseQualityResults mResults;
    private boolean mIsValid;

    public BaseQualityRecalibration(
            final SageConfig config, final IndexedFastaSequenceFile refGenome, final String panelBedFile,
            final List<String> tumorIds, final List<String> tumorBams)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mTumorBams = tumorBams;
        mTumorIds = tumorIds;
        mPanelBedFile = panelBedFile;

        mSampleRecalibrationMap = Maps.newHashMap();
        mRegions = new ConcurrentLinkedQueue<>();
        mResults = new BaseQualityResults();
        mIsValid = true;
    }

    public boolean isValid(){ return mIsValid; }

    public Map<String, BqrRecordMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    public void produceRecalibrationMap()
    {
        if(!mConfig.QualityRecalibration.Enabled)
        {
            buildEmptyRecalibrations();
            return;
        }

        if(mConfig.QualityRecalibration.LoadBqrFiles)
        {
            Map<String,String> sampleFileNames = Maps.newHashMap();
            String outputDir = mConfig.outputDir();
            mConfig.ReferenceIds.forEach(x -> sampleFileNames.put(x, generateBqrFilename(outputDir, x)));
            mTumorIds.forEach(x -> sampleFileNames.put(x, generateBqrFilename(outputDir, x)));

            for(Map.Entry<String,String> entry : sampleFileNames.entrySet())
            {
                String sampleId = entry.getKey();
                String filename = entry.getValue();

                final List<BqrRecord> counts = BqrFile.read(filename);

                if(counts == null)
                {
                    mIsValid = false;
                    return;
                }

                SG_LOGGER.info("loaded sample({}) {} base quality recalibration records from {}", sampleId, counts.size(), filename);
                mSampleRecalibrationMap.put(sampleId, new BqrRecordMap(counts));
            }

            return;
        }

        final List<PartitionTask> regions = createRegions();

        for(int i = 0; i < mConfig.ReferenceIds.size(); i++)
        {
            processSample(mConfig.ReferenceIds.get(i), mConfig.ReferenceBams.get(i), regions);
        }

        for(int i = 0; i < mTumorIds.size(); i++)
        {
            processSample(mTumorIds.get(i), mTumorBams.get(i), regions);
        }

        if(mConfig.logPerfStats())
            mResults.logPerfStats();

        SG_LOGGER.info("base quality recalibration cache generated");
    }

    private void processSample(final String sampleId, final String bamFile, final List<PartitionTask> regions)
    {
        mRegions.addAll(regions);
        mResults.clear();

        int regionCount = mRegions.size();

        SG_LOGGER.debug("samples({}) building base-qual recalibration map from {} regions", sampleId, regionCount);

        List<BqrThread> workers = new ArrayList<>();

        for(int i = 0; i < min(mRegions.size(), mConfig.Threads); ++i)
        {
            workers.add(new BqrThread(mConfig, mRefGenome, bamFile, mRegions, mResults));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                SG_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        // merge results for this sample across all regions
        final Map<BqrKey,Integer> allQualityCounts = mResults.getCombinedQualityCounts();

        final List<BqrRecord> records = convertToRecords(allQualityCounts);

        mSampleRecalibrationMap.put(sampleId, new BqrRecordMap(records));

        // write results to file
        if(mConfig.QualityRecalibration.WriteFile)
            writeSampleData(sampleId, records);
    }

    private void buildEmptyRecalibrations()
    {
        for(String sample : mConfig.ReferenceIds)
        {
            mSampleRecalibrationMap.put(sample, new BqrRecordMap(Collections.emptyList()));
        }

        for(String sample : mTumorIds)
        {
            mSampleRecalibrationMap.put(sample, new BqrRecordMap(Collections.emptyList()));
        }
    }

    private List<BqrRecord> convertToRecords(final Map<BqrKey,Integer> allQualityCounts)
    {
        final List<BqrRecord> result = Lists.newArrayList();

        final Map<BqrKey,Integer> refCountMap = allQualityCounts.entrySet().stream()
                .filter(x -> x.getKey().Ref == x.getKey().Alt)
                .collect(Collectors.toMap(x -> x.getKey(), x -> x.getValue()));

        for(Map.Entry<BqrKey,Integer> entry : allQualityCounts.entrySet())
        {
            BqrKey key = entry.getKey();

            if(key.Quality == 0)
                continue;

            BqrKey refKey = new BqrKey(key.Ref, key.Ref, key.TrinucleotideContext, key.Quality);

            int refCount = refCountMap.getOrDefault(refKey, 0);

            if(refCount > 0)
            {
                double recalibratedQual = key.Alt == key.Ref ? key.Quality : recalibratedQual(refCount, entry.getValue());
                result.add(new BqrRecord(key, entry.getValue(), recalibratedQual));
            }
        }

        return result;
    }

    public static double recalibratedQual(int refCount, int altCount)
    {
        double percent = altCount / (double) (altCount + refCount);
        return -10 * Math.log10(percent);
    }

    private static final int END_BUFFER = 1000000;
    private static final int REGION_SIZE = 100000;

    private List<PartitionTask> createRegions()
    {
        List<PartitionTask> regionTasks = Lists.newArrayList();
        int taskId = 1;

        // form regions from 2MB per chromosome and additionally include the coding panel
        Map<Chromosome,List<BaseRegion>> panelBed = !mPanelBedFile.isEmpty() ? loadBedFileChrMap(mPanelBedFile) : null;

        for(final SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            final String chromosome = sequenceRecord.getSequenceName();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) || !HumanChromosome.fromString(chromosome).isAutosome())
                continue;

            List<PartitionTask> chrRegionTasks = Lists.newArrayList();

            int start = sequenceRecord.getSequenceLength() - END_BUFFER - mConfig.QualityRecalibration.SampleSize;
            int end = sequenceRecord.getSequenceLength() - (END_BUFFER + 1);

            while(start < end)
            {
                chrRegionTasks.add(new PartitionTask(new ChrBaseRegion(chromosome, start, start + REGION_SIZE - 1), taskId++));
                start += REGION_SIZE;
            }

            if(panelBed != null && panelBed.containsKey(HumanChromosome.fromString(chromosome)))
            {
                List<BaseRegion> panelRegions = panelBed.get(HumanChromosome.fromString(chromosome));

                List<PartitionTask> panelRegionTasks = Lists.newArrayList();
                for(BaseRegion panelRegion : panelRegions)
                {
                    if(chrRegionTasks.stream().anyMatch(x -> panelRegion.overlaps(x.Partition)))
                        continue;

                    panelRegionTasks.add(new PartitionTask(new ChrBaseRegion(chromosome, panelRegion.start(), panelRegion.end()), taskId++));
                }

                chrRegionTasks.addAll(panelRegionTasks);
            }

            regionTasks.addAll(chrRegionTasks);
        }

        Collections.sort(regionTasks, new PartitionTask.PartitionTaskComparator());

        return regionTasks;
    }

    private void writeSampleData(final String sampleId, final Collection<BqrRecord> records)
    {
        try
        {
            final String tsvFile = generateBqrFilename(mConfig.outputDir(), sampleId);
            SG_LOGGER.debug("writing base quality recalibration file: {}", tsvFile);


            BqrFile.write(tsvFile, records.stream().collect(Collectors.toList()));

            if(mConfig.QualityRecalibration.WritePlot)
            {
                RExecutor.executeFromClasspath("r/baseQualityRecalibrationPlot.R", tsvFile);
            }
        }
        catch(Exception e)
        {
            SG_LOGGER.error(" sample({}) failed to write base recalibration: {}", sampleId, e.toString());
        }
    }
}
