package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.quality.QualityRecalibrationFile.generateBqrFilename;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BaseQualityRecalibration
{
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;

    private final Map<String,QualityRecalibrationMap> mSampleRecalibrationMap;
    private final Queue<PartitionTask> mRegions;
    private final BaseQualityResults mResults;

    public BaseQualityRecalibration(final SageConfig config, final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
        mConfig = config;

        mSampleRecalibrationMap = Maps.newHashMap();
        mRegions = new ConcurrentLinkedQueue<>();
        mResults = new BaseQualityResults();
     }

    public Map<String,QualityRecalibrationMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    public static Map<String,QualityRecalibrationMap> buildQualityRecalibrationMap(
            final SageConfig config, final IndexedFastaSequenceFile refGenome)
    {
        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(config, refGenome);
        baseQualityRecalibration.produceRecalibrationMap();
        return baseQualityRecalibration.getSampleRecalibrationMap();
    }

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
            String outputDir = mConfig.formOutputDir();
            mConfig.ReferenceIds.forEach(x -> sampleFileNames.put(x, generateBqrFilename(x, outputDir)));
            mConfig.TumorIds.forEach(x -> sampleFileNames.put(x, generateBqrFilename(x, outputDir)));

            for(Map.Entry<String,String> entry : sampleFileNames.entrySet())
            {
                String sampleId = entry.getKey();
                String filename = entry.getValue();

                final List<QualityRecalibrationRecord> counts = QualityRecalibrationFile.read(filename);

                SG_LOGGER.info("loaded sample({}) {} base quality recalibration records from {}", sampleId, counts.size(), filename);
                mSampleRecalibrationMap.put(sampleId, new QualityRecalibrationMap(counts));
            }

            return;
        }

        final List<PartitionTask> regions = createRegions();

        for(int i = 0; i < mConfig.ReferenceIds.size(); i++)
        {
            processSample(mConfig.ReferenceIds.get(i), mConfig.ReferenceBams.get(i), regions);
        }

        for(int i = 0; i < mConfig.TumorIds.size(); i++)
        {
            processSample(mConfig.TumorIds.get(i), mConfig.TumorBams.get(i), regions);
        }

        if(SG_LOGGER.isDebugEnabled())
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
            }
        }

        // merge results for this sample across all regions
        final Map<BaseQualityKey,Integer> allQualityCounts = mResults.getCombinedQualityCounts();

        final List<QualityRecalibrationRecord> records = convertToRecords(allQualityCounts);

        mSampleRecalibrationMap.put(sampleId, new QualityRecalibrationMap(records));

        // write results to file
        if(mConfig.QualityRecalibration.WriteFile)
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

    private List<PartitionTask> createRegions()
    {
        List<PartitionTask> regionTasks = Lists.newArrayList();
        int taskId = 1;

        for(final SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            final String chromosome = sequenceRecord.getSequenceName();

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) || !HumanChromosome.fromString(chromosome).isAutosome())
                continue;

            int start = sequenceRecord.getSequenceLength() - END_BUFFER - mConfig.QualityRecalibration.SampleSize;
            int end = sequenceRecord.getSequenceLength() - (END_BUFFER + 1);

            while(start < end)
            {
                regionTasks.add(new PartitionTask(new ChrBaseRegion(chromosome, start, start + REGION_SIZE - 1), taskId++));
                start += REGION_SIZE;
            }
        }
        return regionTasks;
    }

    private void writeSampleData(final String sampleId, final Collection<QualityRecalibrationRecord> records)
    {
        try
        {
            final String tsvFile = generateBqrFilename(sampleId, mConfig.formOutputDir());
            SG_LOGGER.debug("writing base quality recalibration file: {}", tsvFile);

            QualityRecalibrationFile.write(tsvFile, records.stream().collect(Collectors.toList()));

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
