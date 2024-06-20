package com.hartwig.hmftools.sage.bqr;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.genome.bed.BedFileReader.loadBedFileChrMap;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.common.sage.SageCommon.generateBqrFilename;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.BQR_SAMPLE_SIZE;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.qual.BqrFile;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrRecord;
import com.hartwig.hmftools.common.region.PartitionUtils;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class BaseQualityRecalibration
{
    private final SageConfig mConfig;
    private final String mPanelBedFile;
    private final List<String> mTumorIds;
    private final List<String> mTumorBams;
    private final IndexedFastaSequenceFile mRefGenome;

    private final Map<String,BqrRecordMap> mSampleRecalibrationMap;
    private final Map<String,List<Integer>> mKnownVariantMap;
    private final Queue<ChrBaseRegion> mRegions;
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
        mKnownVariantMap = Maps.newHashMap();
        mRegions = new ConcurrentLinkedQueue<>();
        mResults = new BaseQualityResults();
        mIsValid = true;
    }

    public boolean isValid(){ return mIsValid; }

    public Map<String,BqrRecordMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    public void setKnownVariants(final List<VariantContext> variants)
    {
        for(VariantContext variant : variants)
        {
            if(VariantType.type(variant) != VariantType.SNP)
                continue;

            String chromosome = variant.getContig();

            List<Integer> positions = mKnownVariantMap.get(chromosome);

            if(positions == null)
            {
                positions = Lists.newArrayList();
                mKnownVariantMap.put(chromosome, positions);
            }

            positions.add(variant.getStart());
        }
    }

    public void produceRecalibrationMap()
    {
        if(!mConfig.BQR.Enabled)
        {
            buildEmptyRecalibrations();
            return;
        }

        if(mConfig.BQR.LoadBqrFiles)
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

        final List<ChrBaseRegion> regions = createRegions();

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

    private void processSample(final String sampleId, final String bamFile, final List<ChrBaseRegion> regions)
    {
        mRegions.addAll(regions);
        mResults.clear();

        int regionCount = mRegions.size();

        SG_LOGGER.debug("samples({}) building base-qual recalibration map from {} regions", sampleId, regionCount);

        BqrRecordWriter recordWriter = new BqrRecordWriter(mConfig, sampleId);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(mRegions.size(), mConfig.Threads); ++i)
        {
            workers.add(new BqrThread(mConfig, mRefGenome, bamFile, mRegions, mResults, recordWriter, mKnownVariantMap));
        }

        if(!runThreadTasks(workers))
            System.exit(1);

        // merge results for this sample across all regions
        final Map<BqrKey,Integer> allQualityCounts = mResults.getCombinedQualityCounts();

        final List<BqrRecord> records = convertToRecords(allQualityCounts);

        mSampleRecalibrationMap.put(sampleId, new BqrRecordMap(records));

        // write results to file
        if(mConfig.BQR.WriteFile)
            writeSampleData(sampleId, records);

        recordWriter.close();
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

    public static List<BqrRecord> convertToRecords(final Map<BqrKey,Integer> allQualityCounts)
    {
        List<BqrRecord> result = Lists.newArrayList();

        // collect the ref counts
        Map<BqrKey,Integer> refCountMap = allQualityCounts.entrySet().stream()
                .filter(x -> x.getKey().Ref == x.getKey().Alt)
                .collect(Collectors.toMap(x -> x.getKey(), x -> x.getValue()));

        // make a map of (per-type) trinuc totals across all entries
        Map<BqrKey,Integer> triNucMap = Maps.newHashMap();
        byte noBaseOrQual = 1;

        for(Map.Entry<BqrKey,Integer> entry : allQualityCounts.entrySet())
        {
            BqrKey key = entry.getKey();
            int count = entry.getValue();

            if(key.Quality == 0)
                continue;

            BqrKey triNucKey = new BqrKey(noBaseOrQual, noBaseOrQual, key.TrinucleotideContext, noBaseOrQual, key.ReadType);

            Integer triNucCount = triNucMap.get(triNucKey);
            triNucMap.put(triNucKey, triNucCount != null ? triNucCount + count : count);
        }

        Set<BqrKey> syntheticAltKeys = Sets.newHashSet();

        for(Map.Entry<BqrKey,Integer> entry : allQualityCounts.entrySet())
        {
            BqrKey key = entry.getKey();

            if(key.Quality == 0)
                continue;

            double recalibratedQual = 0;

            BqrKey refKey = new BqrKey(key.Ref, key.Ref, key.TrinucleotideContext, key.Quality, key.ReadType);

            if(key.Alt == key.Ref)
            {
                recalibratedQual = key.Quality;
            }
            else
            {
                recalibratedQual = calcRecalibratedQual(key, entry.getValue(), refCountMap, triNucMap);
            }

            result.add(new BqrRecord(key, entry.getValue(), recalibratedQual));

            // add alt entries for any which sample has no results
            for(int i = 0; i < DNA_BASE_BYTES.length; ++i)
            {
                byte alt = DNA_BASE_BYTES[i];

                if(alt == refKey.Ref)
                    continue;

                BqrKey altKey = new BqrKey(key.Ref, alt, key.TrinucleotideContext, key.Quality, key.ReadType);

                if(!allQualityCounts.containsKey(altKey) && !syntheticAltKeys.contains(altKey))
                {
                    syntheticAltKeys.add(altKey);

                    double recalQualMin = calcRecalibratedQual(altKey, 1, refCountMap, triNucMap);
                    double syntheticQual = max(recalQualMin + 10 * log10(2), key.Quality);

                    result.add(new BqrRecord(altKey, 0, syntheticQual));
                }
            }
        }

        return result;
    }

    private static final byte NO_BASE_OR_QUAL = 1; // value is irrelevant, just to complete a map entry

    private static double calcRecalibratedQual(
            final BqrKey key, final int observedCount, final Map<BqrKey,Integer> refCountMap, final Map<BqrKey,Integer> triNucMap)
    {
        // example: to calc qual of T>C errors at a GTA trinucleotide, with qual 30
        // error_rate = (# all GCA sites / # all GTA sites) * (# qual 30 T>Cs GTA sites) / (# qual 30 T>Cs GTA sites + # qual 30 C refs GCA sites)

        // or in general terms:
        // errorRate = (# alt TN sites / # ref TN sites) * (# alts at ref TN sites) / (# alts at ref TN sites + # refs at alt TN sites)

        byte[] altTriNucContext = new byte[] { key.TrinucleotideContext[0], key.Alt, key.TrinucleotideContext[2] };

        BqrKey altKey = new BqrKey(key.Alt, key.Alt, altTriNucContext, key.Quality, key.ReadType);
        int altRefCount = refCountMap.getOrDefault(altKey, 0);

        BqrKey altTriNucKey = new BqrKey(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, altTriNucContext, NO_BASE_OR_QUAL, key.ReadType);
        int altTriNucCount = triNucMap.getOrDefault(altTriNucKey, 0);

        BqrKey refTriNucKey = new BqrKey(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, key.TrinucleotideContext, NO_BASE_OR_QUAL, key.ReadType);
        int refTriNucCount = triNucMap.getOrDefault(refTriNucKey, 0);

        double triNucRate = refTriNucCount > 0 ? altTriNucCount / (double)refTriNucCount : 0;

        if(triNucRate == 0 || (observedCount + altRefCount) == 0)
            return 0;

        double calcProbability = triNucRate * observedCount / (observedCount + altRefCount);
        return probabilityToPhredQual(calcProbability);
    }

    private static final int END_BUFFER = 1000000;
    private static final int REGION_SIZE = 100000;

    private List<ChrBaseRegion> createRegions()
    {
        List<ChrBaseRegion> regionTasks = Lists.newArrayList();

        if(!mConfig.BQR.UsePanel && !mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
            {
                regionTasks.add(new ChrBaseRegion(
                        region.Chromosome, region.start() - REGION_SIZE, region.end() + REGION_SIZE - 1));
            }

            return regionTasks;
        }

        if(mConfig.BQR.FullBam)
        {
            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                if(chromosome.isAutosome())
                {
                    regionTasks.addAll(PartitionUtils.partitionChromosome(chromosome.toString(), mConfig.RefGenVersion, BQR_SAMPLE_SIZE));
                }
            }

            return regionTasks;
        }

        // form regions from 2MB per chromosome and additionally include the coding panel
        Map<Chromosome,List<BaseRegion>> panelBed = !mPanelBedFile.isEmpty() ? loadBedFileChrMap(mPanelBedFile) : null;

        for(final SAMSequenceRecord sequenceRecord : mRefGenome.getSequenceDictionary().getSequences())
        {
            final String chromosome = sequenceRecord.getSequenceName();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(!HumanChromosome.contains(chromosome) || HumanChromosome.fromString(chromosome).isAllosome())
                continue;

            List<ChrBaseRegion> chrRegionTasks = Lists.newArrayList();

            int start = sequenceRecord.getSequenceLength() - END_BUFFER - mConfig.BQR.SampleSize;
            int end = sequenceRecord.getSequenceLength() - (END_BUFFER + 1);

            while(start < end)
            {
                chrRegionTasks.add(new ChrBaseRegion(chromosome, start, start + REGION_SIZE - 1));
                start += REGION_SIZE;
            }

            if(panelBed != null && panelBed.containsKey(HumanChromosome.fromString(chromosome)))
            {
                List<BaseRegion> panelRegions = panelBed.get(HumanChromosome.fromString(chromosome));

                List<ChrBaseRegion> panelRegionTasks = Lists.newArrayList();
                for(BaseRegion panelRegion : panelRegions)
                {
                    if(chrRegionTasks.stream().anyMatch(x -> panelRegion.overlaps(x)))
                        continue;

                    panelRegionTasks.add(new ChrBaseRegion(chromosome, panelRegion.start(), panelRegion.end()));
                }

                chrRegionTasks.addAll(panelRegionTasks);
            }

            regionTasks.addAll(chrRegionTasks);
        }

        Collections.sort(regionTasks);

        return regionTasks;
    }



    private void writeSampleData(final String sampleId, final Collection<BqrRecord> records)
    {
        try
        {
            final String tsvFile = generateBqrFilename(mConfig.outputDir(), sampleId);
            SG_LOGGER.debug("writing base quality recalibration file: {}", tsvFile);


            BqrFile.write(tsvFile, records.stream().collect(Collectors.toList()));

            if(mConfig.BQR.WritePlot)
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
