package com.hartwig.hmftools.bamtools.utils;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMRecord;

public class BamUmiAnalyser
{
    private final UmiAnalyserConfig mConfig;

    public BamUmiAnalyser(final ConfigBuilder configBuilder)
    {
        mConfig = new UmiAnalyserConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        BT_LOGGER.info("running UMI analyser");

        List<String> knownUmis = loadDelimitedIdFile(mConfig.KnownUmiFile, "KnownUmi", CSV_DELIM);
        Set<String> knownUmiSet = Sets.newHashSet();
        knownUmis.forEach(x -> knownUmiSet.add(x));

        BufferedWriter writer = initialiseWriter();

        List<ChrBaseRegion> partitionRegions = PartitionTask.splitRegionsIntoPartitions(
                mConfig.RawBamFile, mConfig.RefGenomeFile, mConfig.Threads, mConfig.SpecificChrRegions, mConfig.PartitionSize);

        List<PartitionUmiAnalyser> partitionUmiAnalysers = Lists.newArrayList();
        List<Thread> threadTasks = Lists.newArrayList();

        Queue<ChrBaseRegion> regionsQueue = new ConcurrentLinkedQueue<>();
        regionsQueue.addAll(partitionRegions);

        TaskQueue taskQueue = new TaskQueue(regionsQueue, "regions", 100);

        for(int i = 0; i < mConfig.Threads; ++i)
        {
            PartitionUmiAnalyser partitionUmiAnalyser = new PartitionUmiAnalyser(mConfig, taskQueue, knownUmiSet, writer);
            partitionUmiAnalysers.add(partitionUmiAnalyser);
            threadTasks.add(partitionUmiAnalyser);
        }

        BT_LOGGER.debug("splitting {} regions across {} threads", partitionRegions.size(), mConfig.Threads);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        closeBufferedWriter(writer);

        BT_LOGGER.info("UMI analyser complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String filename = mConfig.OutputDir + "umi_analysis.tsv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("ReadId").add("Chromosome").add("PosStart").add("PosEnd").add("Cigar");
            sj.add("MateChr").add("MatePosStart").add("Flags").add("FirstInPair").add("Reversed").add("MateReversed");

            sj.add("UmiReadTag").add("UmiSequence").add("SeqMatched");
            sj.add("RevSequence").add("RevSeqMatched").add("OtherStrippedLength");

            sj.add("BaseLength").add("RawBases");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to initialise writer: {}", e.toString());
            return null;
        }
    }

    protected synchronized static void writeReadInfo(
            final BufferedWriter writer, final SAMRecord rawRead, final String umiReadTag, final String umiSequence,
            boolean seqMatched, final String revUmiSequence, boolean revSeqMatched, int otherStrippedLength)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(rawRead.getReadName());
            sj.add(rawRead.getContig());
            sj.add(String.valueOf(rawRead.getAlignmentStart()));
            sj.add(String.valueOf(rawRead.getAlignmentEnd()));
            sj.add(rawRead.getCigarString());
            sj.add(rawRead.getMateReferenceName());
            sj.add(String.valueOf(rawRead.getMateAlignmentStart()));
            sj.add(String.valueOf(rawRead.getFlags()));

            sj.add(String.valueOf(rawRead.getFirstOfPairFlag()));
            sj.add(String.valueOf(rawRead.getReadNegativeStrandFlag()));
            sj.add(String.valueOf(rawRead.getMateNegativeStrandFlag()));

            sj.add(umiReadTag);
            sj.add(umiSequence);
            sj.add(String.valueOf(seqMatched));

            sj.add(revUmiSequence);
            sj.add(String.valueOf(revSeqMatched));
            sj.add(String.valueOf(otherStrippedLength));

            String readBases = rawRead.getReadString();
            sj.add(String.valueOf(readBases.length()));

            if(rawRead.getReadNegativeStrandFlag())
                readBases = Nucleotides.reverseComplementBases(readBases);

            sj.add(readBases);

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to write read info: {}", e.toString());
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        UmiAnalyserConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamUmiAnalyser bamUmiAnalyser = new BamUmiAnalyser(configBuilder);
        bamUmiAnalyser.run();
    }
}
