package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.FutureTask;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BaseQualityRecalibration
{
    private final ExecutorService mExecutorService;
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;

    private final Map<String,QualityRecalibrationMap> mSampleRecalibrationMap;

    final Map<QualityCounterKey,QualityCounter> mCounters;

    public BaseQualityRecalibration(
            final SageConfig config, final ExecutorService executorService, final IndexedFastaSequenceFile refGenome)
    {
        mExecutorService = executorService;
        mRefGenome = refGenome;
        mConfig = config;

        mSampleRecalibrationMap = Maps.newHashMap();

        mCounters = new ConcurrentHashMap<>(Maps.newHashMap());
    }

    public Map<String,QualityRecalibrationMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    public void produceRecalibrationMap()
    {
        if(!mConfig.QualityRecalibration.Enabled)
        {
            buildEmptyRecalibrations();
            return;
        }

        final List<BaseRegion> regions = createRegions();

        for(int i = 0; i < mConfig.ReferenceIds.size(); i++)
        {
            processSample(mConfig.ReferenceIds.get(i), mConfig.ReferenceBams.get(i), regions);
        }

        for(int i = 0; i < mConfig.TumorIds.size(); i++)
        {
            processSample(mConfig.TumorIds.get(i), mConfig.TumorBams.get(i), regions);
        }
    }

    private void processSample(final String sampleId, final String bamFile, final List<BaseRegion> regions)
    {
        List<BaseQualityRegionCounter> regionCounters = Lists.newArrayList();
        List<FutureTask> taskList = new ArrayList<FutureTask>();

        for(BaseRegion region : regions)
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

        // merge results for this sample across all regions
        final Map<QualityCounterKey,QualityCounter> allQualityCounts = Maps.newHashMap();

        for(BaseQualityRegionCounter regionCounter : regionCounters)
        {
            mergeQualityCounts(allQualityCounts, regionCounter.getQualityCounts());
        }

        final List<QualityRecalibrationRecord> records = convertToRecords(allQualityCounts);

        mSampleRecalibrationMap.put(sampleId, new QualityRecalibrationMap(records));

        // write results to file
        writeSampleData(sampleId, records);
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

    private void mergeQualityCounts(final Map<QualityCounterKey,QualityCounter> allQualityCounts, final Collection<QualityCounter> qualityCounts)
    {
        for(QualityCounter count : qualityCounts)
        {
            QualityCounterKey key = withoutPosition(count);
            QualityCounter counter = allQualityCounts.get(key);
            if(counter == null)
            {
                counter = new QualityCounter(key);
                allQualityCounts.put(key, counter);
            }

            counter.increment(count.count());
        }
    }

    private List<QualityRecalibrationRecord> convertToRecords(final Map<QualityCounterKey,QualityCounter> allQualityCounts)
    {
        final List<QualityCounter> sortedList = Lists.newArrayList(allQualityCounts.values());
        Collections.sort(sortedList);

        final List<QualityRecalibrationRecord> result = Lists.newArrayList();

        final Map<QualityRecalibrationKey,Integer> refCountMap = sortedList.stream()
                .filter(x -> x.ref() == x.alt())
                .collect(Collectors.toMap(x -> QualityRecalibrationKey.from(x.Key), x -> x.count()));

        for(QualityCounter cleanedRecord : sortedList)
        {
            final QualityRecalibrationKey key = QualityRecalibrationKey.from(cleanedRecord.Key);
            final QualityRecalibrationKey refKey = new QualityRecalibrationKey(key.Ref, cleanedRecord.ref(), key.TrinucleotideContext, key.Quality);

            int refCount = refCountMap.getOrDefault(refKey, 0);
            if(refCount > 0)
            {
                double recalibratedQual = cleanedRecord.alt() == cleanedRecord.ref()
                        ? cleanedRecord.qual()
                        : recalibratedQual(refCount, cleanedRecord.count());

                result.add(new QualityRecalibrationRecord(key, cleanedRecord.count(), recalibratedQual));
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

    private List<BaseRegion> createRegions()
    {
        List<BaseRegion> result = Lists.newArrayList();

        // mRefGenome, mConfig.Chromosomes, mConfig.QualityRecalibration.SampleSize
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
                result.add(new BaseRegion(chromosome, start, start + REGION_SIZE - 1));
                start += REGION_SIZE;
            }
        }
        return result;
    }

    private static QualityCounterKey withoutPosition(final QualityCounter count)
    {
        return new QualityCounterKey(count.ref(), count.alt(), count.qual(), 0, count.trinucleotideContext());
    }

    private void writeSampleData(final String sampleId, final Collection<QualityRecalibrationRecord> records)
    {
        try
        {
            final String tsvFile = mConfig.baseQualityRecalibrationFile(sampleId);
            SG_LOGGER.debug("writing base quality recalibration file: {}", tsvFile);

            QualityRecalibrationFile.write(tsvFile, records);

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
