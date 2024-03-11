package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.APP_NAME;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FragmentMetrics
{
    private final MetricsConfig mConfig;

    public FragmentMetrics(final ConfigBuilder configBuilder)
    {
        mConfig = new MetricsConfig(configBuilder);
    }

    public void run()
    {
        CB_LOGGER.info("sample({}) starting fragment GC metrics", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        List<ChrBaseRegion> allRegions = Lists.newArrayList();

        if(!mConfig.TargetRegions.isEmpty())
        {
            int targetRegionCount = mConfig.TargetRegions.values().stream().mapToInt(x -> x.size()).sum();

            int targetRegionBaseCount = mConfig.TargetRegions.values().stream()
                    .mapToInt(x -> x.stream().mapToInt(y -> y.baseLength()).sum()).sum();

            CB_LOGGER.info("capturing data for {} target regions, total base count({})",
                    targetRegionCount, targetRegionBaseCount);
        }

        if(!mConfig.OnlyTargetRegions)
        {
            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                    continue;

                allRegions.addAll(partitionChromosome(
                        chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize));
            }
        }
        else
        {
            // consider expanding partitions to cover regions within range of each other
            for(Map.Entry<String,List<BaseRegion>> entry : mConfig.TargetRegions.entrySet())
            {
                entry.getValue().forEach(x -> allRegions.add(new ChrBaseRegion(entry.getKey(), x.start(), x.end())));
            }

            Collections.sort(allRegions);
        }

        List<TargetRegionData> targetRegionData = Lists.newArrayList();

        Queue<Partition> partitions = new ConcurrentLinkedQueue<>();

        for(int i = 0; i < allRegions.size(); ++i)
        {
            ChrBaseRegion region = allRegions.get(i);

            Partition partition = new Partition(region);

            List<BaseRegion> targetRegions = mConfig.TargetRegions.get(region.Chromosome);

            if(targetRegions != null)
            {
                // if a region spans 2 partition then it will be recorded by each one and results written twice
                List<TargetRegionData> overlappingRegions = targetRegions.stream().filter(x -> x.overlaps(region))
                        .map(x -> new TargetRegionData(region.Chromosome, x)).collect(Collectors.toList());

                partition.TargetRegions.addAll(overlappingRegions);

                targetRegionData.addAll(overlappingRegions);
            }

            partitions.add(partition);
        }

        List<PartitionReader> partitionReaders = Lists.newArrayList();

        if(allRegions.size() == 1 || mConfig.Threads <= 1)
        {
            PartitionReader partitionReader = new PartitionReader(mConfig, partitions);
            partitionReader.run();
            partitionReaders.add(partitionReader);
        }
        else
        {
            CB_LOGGER.info("splitting {} regions across {} threads", allRegions.size(), mConfig.Threads);

            List<Thread> workers = new ArrayList<>();

            for(int i = 0; i < min(allRegions.size(), mConfig.Threads); ++i)
            {
                PartitionReader partitionReader = new PartitionReader(mConfig, partitions);
                partitionReader.start();

                partitionReaders.add(partitionReader);
                workers.add(partitionReader);
            }

            if(!runThreadTasks(workers))
                System.exit(1);
        }

        CB_LOGGER.debug("all regions complete");

        FragmentGcMap combinedAllFragGcMap = new FragmentGcMap();
        FragmentGcMap combinedTargedFragGcMap = new FragmentGcMap();
        FragmentGcMap combinedNonTargedFragGcMap = new FragmentGcMap();

        for(PartitionReader partitionReader : partitionReaders)
        {
            combinedAllFragGcMap.merge(partitionReader.allFragmentGcMap());
            combinedTargedFragGcMap.merge(partitionReader.targetedFragmentGcMap());
            combinedNonTargedFragGcMap.merge(partitionReader.nonTargetedFragmentGcMap());
        }

        writeResults(combinedAllFragGcMap, combinedTargedFragGcMap, combinedNonTargedFragGcMap, targetRegionData);

        CB_LOGGER.info("FragmentGcMetrics complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void writeResults(
            final FragmentGcMap combinedAllFragGcMap, final FragmentGcMap combinedTargedFragGcMap,
            final FragmentGcMap combinedNonTargedFragGcMap, final List<TargetRegionData> targetRegionData)
    {
        try
        {
            // write summary metrics
            String filename = mConfig.OutputDir + mConfig.SampleId + ".fragment_gc_metrics";

            if(mConfig.OutputId != null)
                filename += "." + mConfig.OutputId;

            filename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sb = new StringJoiner(TSV_DELIM);

            if(mConfig.CaptureRegionCounts)
            {
                sb.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END);
            }
            else
            {
                sb.add("RegionType");
            }

            sb.add("FragmentLength").add("DuplicateCount").add("GcPercent").add("Count");
            writer.write(sb.toString());
            writer.newLine();

            combinedAllFragGcMap.writeMap(writer, "ALL");
            combinedTargedFragGcMap.writeMap(writer, "TARGET");
            combinedNonTargedFragGcMap.writeMap(writer, "NON_TARGET");

            if(mConfig.CaptureRegionCounts)
            {
                for(TargetRegionData targetRegion : targetRegionData)
                {
                    String regionStr = format("%s\t%d\t%d", targetRegion.Chromosome, targetRegion.start(), targetRegion.end());
                    targetRegion.FragmentGcCounts.writeMap(writer, regionStr);
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write fragment GC counts file: {}", e.toString());
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        MetricsConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FragmentMetrics fragmentMetrics = new FragmentMetrics(configBuilder);
        fragmentMetrics.run();
    }
}
