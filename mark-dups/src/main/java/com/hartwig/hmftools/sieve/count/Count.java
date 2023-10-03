package com.hartwig.hmftools.sieve.count;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class Count
{
    private final CountConfig mConfig;

    public Count(final ConfigBuilder configBuilder)
    {
        mConfig = new CountConfig(configBuilder);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        CountConfig.addConfig(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        // TODO(m_cooper): Fill in.
        // logVersion();

        Count count = new Count(configBuilder);
        count.run();
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            MD_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        if(mConfig.BucketSize <= 0)
        {
            MD_LOGGER.error("bucket_size must be a positive integer.");
            System.exit(1);
        }

        final Set<String> humanChromosomeSet =
                Arrays.stream(HumanChromosome.values()).map(chr -> stripChrPrefix(chr.toString())).collect(Collectors.toSet());
        final RefGenomeSource refGenome = loadRefGenome(mConfig.RefGenome);
        SAMSequenceDictionary seqDict = refGenome.mRefGenome.getSequenceDictionary();
        if(seqDict == null)
        {
            MD_LOGGER.error("Cannot find .dict file for ref genome {}", mConfig.RefGenome);
            System.exit(1);
        }

        final Vector<CountResult> results = new Vector<>();
        final List<CountTask> countTasks = new ArrayList<>();
        for(SAMSequenceRecord seq : seqDict.getSequences())
        {
            CountTask countTask = null;
            final String contig = seq.getContig();
            if(humanChromosomeSet.contains(stripChrPrefix(contig)))
            {
                countTask = new CountTask(mConfig, true, seq, results);
            }
            else if(contig.endsWith("_alt"))
            {
                countTask = new CountTask(mConfig, false, seq, results);
            }
            else if(contig.endsWith("_decoy"))
            {
                countTask = new CountTask(mConfig, false, seq, results);
            }
            else if(contig.equals("chrEBV"))
            {
                countTask = new CountTask(mConfig, false, seq, results);
            }

            if(countTask != null)
            {
                countTasks.add(countTask);
            }
        }

        final List<Callable> callableList = countTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        writeContigCountOutputFile(results);
        writeBucketCountOutputFile(results);

        MD_LOGGER.info("count complete");
    }

    private void writeContigCountOutputFile(final Vector<CountResult> results)
    {
        MD_LOGGER.info("Writing contig count output to {}.", mConfig.ContigCountOutputFile);
        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(mConfig.ContigCountOutputFile));
            writer.write("Contig\tPrimary Reads\tSupplementary Reads");
            writer.newLine();
            for(CountResult result : results)
            {
                writer.write(result.getSeq().getContig() + '\t' + result.getPrimaryReadCount() + '\t' + result.getSuppReadCount());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("An exception was raised while writing contig output to {}: {}", mConfig.ContigCountOutputFile, e.toString());
            System.exit(1);
        }
    }

    private void writeBucketCountOutputFile(final Vector<CountResult> results)
    {
        MD_LOGGER.info("Writing bucket count output to {}.", mConfig.BucketCountOutputFile);

        final TreeSet<Integer> allMapQ = new TreeSet<>();
        for(CountResult result : results)
        {
            List<Map<Integer, Long>> buckets = result.getPrimaryBucketCounts();
            if(buckets != null)
            {
                for(Map<Integer, Long> bucket : buckets)
                {
                    allMapQ.addAll(bucket.keySet());
                }
            }

            buckets = result.getSuppBucketCounts();
            if(buckets != null)
            {
                for(Map<Integer, Long> bucket : buckets)
                {
                    allMapQ.addAll(bucket.keySet());
                }
            }
        }

        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(mConfig.BucketCountOutputFile));
            writer.write("Contig\tStartPos\tEndPos");
            for(int mapQ : allMapQ)
            {
                writer.write("\tPrimary MAPQ " + mapQ);
            }

            for(int mapQ : allMapQ)
            {
                writer.write("\tSupp MAPQ " + mapQ);
            }

            writer.newLine();

            for(CountResult result : results)
            {
                final String contig = result.getSeq().getContig();
                final int seqLength = result.getSeq().getSequenceLength();
                final List<Map<Integer, Long>> primaryBuckets = result.getPrimaryBucketCounts();
                if(primaryBuckets == null)
                {
                    continue;
                }

                final List<Map<Integer, Long>> suppBuckets = result.getSuppBucketCounts();
                if(suppBuckets == null)
                {
                    continue;
                }

                for(int i = 0; i < primaryBuckets.size(); i++)
                {
                    final int startPos = i * mConfig.BucketSize + 1;
                    final int endPos = Math.min((i + 1) * mConfig.BucketSize, seqLength);
                    writer.write(contig + "\t" + startPos + "\t" + endPos);
                    for(int mapQ : allMapQ)
                    {
                        writer.write("\t" + primaryBuckets.get(i).getOrDefault(mapQ, (long) 0));
                    }

                    for(int mapQ : allMapQ)
                    {
                        writer.write("\t" + suppBuckets.get(i).getOrDefault(mapQ, (long) 0));
                    }

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            MD_LOGGER.error("An exception was raised while writing bucket count output to {}: {}", mConfig.BucketCountOutputFile, e.toString());
            System.exit(1);
        }
    }
}
