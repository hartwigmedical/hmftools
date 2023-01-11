package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class MarkDuplicates
{
    private final MarkDupsConfig mConfig;

    public MarkDuplicates(final CommandLine cmd)
    {
        mConfig = new MarkDupsConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BM_LOGGER.info("sample({}) starting mark duplicates", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        List<ChromosomeReader> chromosomeReaders = Lists.newArrayList();

        RefGenomeCoordinates refGenomeCoordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        RecordWriter recordWriter = new RecordWriter(mConfig);
        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);
        final List<Callable> callableList = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            ChrBaseRegion chrBaseRegion = new ChrBaseRegion(chromosomeStr, 1, refGenomeCoordinates.Lengths.get(chromosome));

            ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, recordWriter, partitionDataStore);
            chromosomeReaders.add(chromosomeReader);
            callableList.add(chromosomeReader);
        }

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            List<Fragment> fragments = partitionData.extractRemainingFragments();
            recordWriter.writeFragments(fragments);
        }

        if(mConfig.UseInterimFiles)
        {
            recordWriter.closeInterimFiles();

            // routine for handling these..
        }

        recordWriter.close();

        long totalProcessReads = chromosomeReaders.stream().mapToLong(x -> x.totalRecordCount()).sum();
        BM_LOGGER.debug("all chromosome tasks complete, reads processed({})", totalProcessReads);

        DuplicateStats combinedStats = new DuplicateStats();
        chromosomeReaders.forEach(x -> combinedStats.merge(x.duplicateStats()));

        BM_LOGGER.info("stats: duplicates({})", combinedStats.Duplicates);

        if(BM_LOGGER.isDebugEnabled())
        {
            List<Integer> frequencies = combinedStats.DuplicateFrequencies.keySet().stream().collect(Collectors.toList());
            Collections.sort(frequencies);

            for(Integer frequency : frequencies)
            {
                BM_LOGGER.debug("duplicate frequency({}={})", frequency, combinedStats.DuplicateFrequencies.get(frequency));
            }
        }

        if(totalProcessReads != recordWriter.recordWriteCount())
        {
            BM_LOGGER.warn("reads processed({}) vs written({}) mismatch", totalProcessReads, recordWriter.recordWriteCount());
            checkMissingReads(chromosomeReaders, recordWriter);
        }

        PerformanceCounter combinedPerfCounter = chromosomeReaders.get(0).perfCounter();

        for(int i = 1; i < chromosomeReaders.size(); ++i)
        {
            combinedPerfCounter.merge(chromosomeReaders.get(i).perfCounter());
        }

        if(mConfig.PerfDebug)
            combinedPerfCounter.logIntervalStats(10);
        else
            combinedPerfCounter.logStats();

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        BM_LOGGER.info("Mark duplicates complete, mins({})", format("%.3f", timeTakeMins));
    }

    private void checkMissingReads(final List<ChromosomeReader> chromosomeReaders, final RecordWriter recordWriter)
    {
        if(!mConfig.runReadChecks())
            return;

        Set<SAMRecord> recordsProcessed = Sets.newHashSet();
        chromosomeReaders.forEach(x -> recordsProcessed.addAll(x.readsProcessed()));

        List<SAMRecord> recordsWritten = recordWriter.readsWritten().stream().collect(Collectors.toList());

        for(SAMRecord readProcessed : recordsProcessed)
        {
            boolean found = false;
            for(int i = 0; i < recordsWritten.size(); ++i)
            {
                SAMRecord read = recordsWritten.get(i);

                if(read == readProcessed)
                {
                    recordsWritten.remove(i);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                BM_LOGGER.error("read processed but not written: {}", readToString(readProcessed));
            }
        }

        if(!recordsWritten.isEmpty())
        {
            for(int i = 0; i < recordsWritten.size(); ++i)
            {
                SAMRecord read = recordsWritten.get(i);

                BM_LOGGER.error("read extra written: {}", readToString(read));
            }
        }
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("bam-tools.version");
        BM_LOGGER.info("BamTools version: {}", version.version());

        final Options options = MarkDupsConfig.createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            MarkDuplicates markDuplicates = new MarkDuplicates(cmd);
            markDuplicates.run();
        }
        catch(ParseException e)
        {
            BM_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("MarkDuplicates", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
