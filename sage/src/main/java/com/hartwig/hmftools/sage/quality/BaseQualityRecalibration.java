package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.FutureTask;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BaseQualityRecalibration
{
    private final ExecutorService mExecutorService;
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;

    private final Map<String,QualityRecalibrationMap> mSampleRecalibrationMap;

    private final PerformanceCounter mPerfCounter;

    public BaseQualityRecalibration(
            final SageConfig config, final ExecutorService executorService, final IndexedFastaSequenceFile refGenome)
    {
        mExecutorService = executorService;
        mRefGenome = refGenome;
        mConfig = config;

        mSampleRecalibrationMap = Maps.newHashMap();

        mPerfCounter = new PerformanceCounter("BaseQualRecal");
     }

    public Map<String,QualityRecalibrationMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    public void produceRecalibrationMap()
    {
        if(!mConfig.QualityRecalibration.Enabled)
        {
            buildEmptyRecalibrations();
            return;
        }

        final List<ChrBaseRegion> regions = createRegions();

        for(int i = 0; i < mConfig.ReferenceIds.size(); i++)
        {
            processSample(mConfig.ReferenceIds.get(i), mConfig.ReferenceBams.get(i), regions);
        }

        if(SG_LOGGER.isDebugEnabled())
            mPerfCounter.logStats();

        for(int i = 0; i < mConfig.TumorIds.size(); i++)
        {
            processSample(mConfig.TumorIds.get(i), mConfig.TumorBams.get(i), regions);
        }

        if(SG_LOGGER.isDebugEnabled())
            mPerfCounter.logStats();
    }

    private void processSample(final String sampleId, final String bamFile, final List<ChrBaseRegion> regions)
    {
        List<BaseQualityRegionCounter> regionCounters = Lists.newArrayList();
        List<FutureTask> taskList = new ArrayList<FutureTask>();

        for(ChrBaseRegion region : regions)
        {
            BaseQualityRegionCounter regionCounter = new BaseQualityRegionCounter(mConfig, bamFile, mRefGenome, region);
            regionCounters.add(regionCounter);

            FutureTask futureTask = new FutureTask(regionCounter);
            taskList.add(futureTask);
            mExecutorService.execute(futureTask);
        }

        try
        {
            for(FutureTask task : taskList)
            {
                task.get();
            }
        }
        catch (Exception e)
        {
            SG_LOGGER.error("sampleId({}) failed to produce base-qual counts: {}", e.toString());
            e.printStackTrace();
        }

        // regionCounters.forEach(x -> x.logPerfs());

        mPerfCounter.start();

        // merge results for this sample across all regions
        final Map<BaseQualityKey,Integer> allQualityCounts = Maps.newHashMap();

        for(BaseQualityRegionCounter regionCounter : regionCounters)
        {
            mergeQualityCounts(allQualityCounts, regionCounter.getQualityCounts());
        }

        final List<QualityRecalibrationRecord> records = convertToRecords(allQualityCounts);

        mSampleRecalibrationMap.put(sampleId, new QualityRecalibrationMap(records));

        // write results to file
        writeSampleData(sampleId, records);

        mPerfCounter.stop();
    }

    private void buildEmptyRecalibrations()
    {
        for(String sample : mConfig.ReferenceIds)
        {
            mSampleRecalibrationMap.put(sample, new QualityRecalibrationMap(Collections.emptyList()));
        }

        for(String sample : mConfig.TumorIds)
        {
            mSampleRecalibrationMap.put(sample, new QualityRecalibrationMap(Collections.emptyList()));
        }
    }

    private void mergeQualityCounts(final Map<BaseQualityKey,Integer> allQualityCounts, final Collection<QualityCounter> qualityCounts)
    {
        for(QualityCounter counter : qualityCounts)
        {
            BaseQualityKey key = counter.Key; //withoutPosition(counter);
            Integer count = allQualityCounts.get(key);

            allQualityCounts.put(key, count != null ? count + counter.count() : counter.count());
        }
    }

    private List<QualityRecalibrationRecord> convertToRecords(final Map<BaseQualityKey,Integer> allQualityCounts)
    {
        final List<QualityRecalibrationRecord> result = Lists.newArrayList();

        final Map<BaseQualityKey,Integer> refCountMap = allQualityCounts.entrySet().stream()
                .filter(x -> x.getKey().Ref == x.getKey().Alt)
                .collect(Collectors.toMap(x -> x.getKey(), x -> x.getValue()));

        for(Map.Entry<BaseQualityKey,Integer> entry : allQualityCounts.entrySet())
        {
            BaseQualityKey key = entry.getKey();
            BaseQualityKey refKey = new BaseQualityKey(key.Ref, key.Ref, key.TrinucleotideContext, key.Quality);

            int refCount = refCountMap.getOrDefault(refKey, 0);
            if(refCount > 0)
            {
                double recalibratedQual = key.Alt == key.Ref
                        ? key.Quality : recalibratedQual(refCount, entry.getValue());

                result.add(new QualityRecalibrationRecord(key, entry.getValue(), recalibratedQual));
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

    private List<ChrBaseRegion> createRegions()
    {
        List<ChrBaseRegion> result = Lists.newArrayList();

        for(final SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            final String chromosome = sequenceRecord.getSequenceName();

            if(!mConfig.Chromosomes.isEmpty() && !mConfig.Chromosomes.contains(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) || !HumanChromosome.fromString(chromosome).isAutosome())
                continue;

            int start = sequenceRecord.getSequenceLength() - END_BUFFER - mConfig.QualityRecalibration.SampleSize;
            int end = sequenceRecord.getSequenceLength() - (END_BUFFER + 1);

            while(start < end)
            {
                result.add(new ChrBaseRegion(chromosome, start, start + REGION_SIZE - 1));
                start += REGION_SIZE;
            }
        }
        return result;
    }

    private void writeSampleData(final String sampleId, final Collection<QualityRecalibrationRecord> records)
    {
        try
        {
            final String tsvFile = mConfig.baseQualityRecalibrationFile(sampleId);
            SG_LOGGER.debug("writing base quality recalibration file: {}", tsvFile);

            QualityRecalibrationFile.write(tsvFile, records.stream().collect(Collectors.toList()));

            if(mConfig.QualityRecalibration.Plot)
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
