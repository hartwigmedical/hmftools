package com.hartwig.hmftools.common.blastn;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class BlastnRunner
{
    // common config
    public static final String BLAST_TOOL = "blast";
    public static final String BLAST_TOOL_DESC = "Path to BlastN";

    public static final String BLAST_DB = "blast_db";
    public static final String BLAST_DB_DESC = "blast_db";

    public static void registerBlastn(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPath(BLAST_TOOL, required, BLAST_TOOL_DESC);
        configBuilder.addPath(BLAST_DB, required, BLAST_DB_DESC);
    }

    private static final Logger sLogger = LogManager.getLogger(BlastnRunner.class);
    private static final String DEFAULT_TASK = "blastn";
    private static final int OUTPUT_STREAM_BUFFER_SIZE = 65536;

    public enum BlastColumns
    {
        qseqid,
        qlen,
        sseqid,
        stitle,
        pident,
        qcovs,
        length,
        mismatch,
        gapopen,
        qstart,
        qend,
        sstart,
        send,
        qframe,
        sframe,
        evalue,
        bitscore,
        qseq,
        sseq
    }

    private final String task;
    private final String prefix;
    private final String blastDir;
    private final String blastDb;
    private final String outputDir;
    private final int numThreads;
    private final boolean keepOutput;
    private final OptionalDouble expectedValueCutoff;
    private final OptionalInt wordSize;
    private final OptionalInt reward;
    private final OptionalInt penalty;
    private final OptionalInt gapopen;
    private final OptionalInt gapextend;
    private final boolean subjectBestHit;

    private BlastnRunner(Builder builder)
    {
        this.task = builder.task;
        this.prefix = builder.prefix;
        this.blastDir = builder.blastDir;
        this.blastDb = builder.blastDb;
        this.outputDir = builder.outputDir;
        this.numThreads = builder.numThreads;
        this.keepOutput = builder.keepOutput;
        this.expectedValueCutoff = builder.expectedValueCutoff;
        this.wordSize = builder.wordSize;
        this.reward = builder.reward;
        this.penalty = builder.penalty;
        this.gapopen = builder.gapopen;
        this.gapextend = builder.gapextend;
        this.subjectBestHit = builder.subjectBestHit;
    }

    @SuppressWarnings("ResultOfMethodCallIgnored")
    public Multimap<Integer, BlastnMatch> run(Map<Integer, String> querySequences)
    {
        if (querySequences == null || querySequences.isEmpty())
        {
            return ArrayListMultimap.create();
        }

        long startTimeMs = System.currentTimeMillis();

        String fastaFile = outputDir + "/" + prefix + ".blastn.fa";
        writeBlastFasta(querySequences, fastaFile);

        List<String> command = Lists.newArrayList(blastDir + "/bin/blastn",
                "-db", "GCF_000001405.39_top_level",
                "-task", task,
                "-mt_mode", "0",
                "-num_threads", Integer.toString(numThreads),
                "-query", fastaFile);

        // Add the optional parameter override values
        expectedValueCutoff.ifPresent(o -> {
            command.add("-evalue");
            command.add(Double.toString(o));
        });
        wordSize.ifPresent(o -> {
            command.add("-word_size");
            command.add(Integer.toString(o));
        });
        reward.ifPresent(o -> {
            command.add("-reward");
            command.add(Integer.toString(o));
        });
        penalty.ifPresent(o -> {
            command.add("-penalty");
            command.add(Integer.toString(o));
        });
        gapopen.ifPresent(o -> {
            command.add("-gapopen");
            command.add(Integer.toString(o));
        });
        gapextend.ifPresent(o -> {
            command.add("-gapextend");
            command.add(Integer.toString(o));
        });
        if(subjectBestHit)
        {
            command.add("-subject_besthit");
        }

        File outputFileCsv = new File(outputDir + "/" + prefix + ".blastn.csv.gz");

        sLogger.debug("running blastn on sample {}, {} sequences, output: {}", prefix, querySequences.size(), outputFileCsv);

        command.add("-outfmt");
        command.add("6 " + Arrays.stream(BlastColumns.values()).map(BlastColumns::name).collect(Collectors.joining(" ")));

        ProcessBuilder processBuilder = new ProcessBuilder(command).redirectError(ProcessBuilder.Redirect.INHERIT);
        Map<String, String> environment = processBuilder.environment();
        environment.put("BLASTDB", blastDb);

        sLogger.debug("{}", String.join(" ", processBuilder.command()));

        try
        {
            Process process = processBuilder.start();
            try(GZIPOutputStream outputStream = new GZIPOutputStream(new FileOutputStream(outputFileCsv)))
            {
                // write the column headers
                String headers = Arrays.stream(BlastColumns.values()).map(BlastColumns::name).collect(Collectors.joining("\t")) + '\n';
                outputStream.write(headers.getBytes(StandardCharsets.UTF_8));

                try(InputStream inputStream = process.getInputStream())
                {
                    byte[] buffer = new byte[OUTPUT_STREAM_BUFFER_SIZE];
                    int read;
                    while ((read = inputStream.read(buffer, 0, OUTPUT_STREAM_BUFFER_SIZE)) >= 0)
                    {
                        outputStream.write(buffer, 0, read);
                    }
                }
            }

            int result = process.waitFor();
            if(result != 0)
            {
                sLogger.fatal("Error executing blastn");
                throw new RuntimeException("blastn execution failed");
            }
        }
        catch(IOException | InterruptedException e)
        {
            sLogger.error("Error running blastn", e);
            throw new RuntimeException(e);
        }

        sLogger.trace("blastn run complete, mins({})", runTimeMinsStr(startTimeMs));

        Multimap<Integer, BlastnMatch> blastnMatches = processBlast(outputFileCsv.getAbsolutePath());

        if(!keepOutput)
        {
            outputFileCsv.delete();
            new File(fastaFile).delete();
        }

        return blastnMatches;
    }

    public static void writeBlastFasta(Map<Integer, String> vdjSequences, String fastaPath)
    {
        try(BufferedWriter writer = FileWriterUtils.createBufferedWriter(fastaPath))
        {
            for(Map.Entry<Integer, String> entry : vdjSequences.entrySet())
            {
                Integer key = entry.getKey();
                String vdjSeq = entry.getValue();
                writer.write(">" + key + "\n");
                writer.write(vdjSeq + "\n");
            }
        }
        catch(IOException e)
        {
            sLogger.error("Error writing FASTA file", e);
            throw new RuntimeException(e);
        }
    }

    public static Multimap<Integer, BlastnMatch> processBlast(String blastOutputCsv)
    {
        Multimap<Integer, BlastnMatch> blastnResults = ArrayListMultimap.create();

        try(DelimFileReader reader = new DelimFileReader(blastOutputCsv))
        {
            for(DelimFileReader.Row row : reader)
            {
                int qframe = row.getInt(BlastColumns.qframe);
                if(qframe != 1)
                {
                    throw new IllegalArgumentException("query frame should always be positive");
                }

                int qseqid = row.getInt(BlastColumns.qseqid);
                BlastnMatch blastnMatch = new BlastnMatch(row.getInt(BlastColumns.qlen),
                        row.get(BlastColumns.stitle),
                        row.getDouble(BlastColumns.pident),
                        row.getDouble(BlastColumns.qcovs),
                        row.getInt(BlastColumns.length),
                        row.getInt(BlastColumns.mismatch),
                        row.getInt(BlastColumns.gapopen),
                        row.getInt(BlastColumns.qstart),
                        row.getInt(BlastColumns.qend),
                        row.getInt(BlastColumns.sstart),
                        row.getInt(BlastColumns.send),
                        Strand.valueOf(row.getInt(BlastColumns.sframe)),
                        row.getDouble(BlastColumns.evalue),
                        row.getDouble(BlastColumns.bitscore),
                        row.get(BlastColumns.qseq),
                        row.get(BlastColumns.sseq));

                blastnResults.put(qseqid, blastnMatch);
            }
        }

        return blastnResults;
    }

    public static class Builder
    {
        private String task = DEFAULT_TASK;
        private String prefix;
        private String blastDir;
        private String blastDb;
        private String outputDir;
        private int numThreads = 1;
        private boolean keepOutput = false;
        private OptionalDouble expectedValueCutoff = OptionalDouble.empty();
        private OptionalInt wordSize = OptionalInt.empty();
        private OptionalInt reward = OptionalInt.empty();
        private OptionalInt penalty = OptionalInt.empty();
        private OptionalInt gapopen = OptionalInt.empty();
        private OptionalInt gapextend = OptionalInt.empty();

        // default subject best hit to true
        private boolean subjectBestHit = true;

        public Builder withTask(String task)
        {
            this.task = task;
            return this;
        }

        public Builder withPrefix(String prefix)
        {
            this.prefix = prefix;
            return this;
        }

        public Builder withBlastDir(String blastDir)
        {
            this.blastDir = blastDir;
            return this;
        }

        public Builder withBlastDb(String blastDb)
        {
            this.blastDb = blastDb;
            return this;
        }

        public Builder withOutputDir(String outputDir)
        {
            this.outputDir = outputDir;
            return this;
        }

        public Builder withNumThreads(int numThreads)
        {
            this.numThreads = numThreads;
            return this;
        }

        public Builder withKeepOutput(boolean keepOutput)
        {
            this.keepOutput = keepOutput;
            return this;
        }

        public Builder withExpectedValueCutoff(double expectedValueCutoff)
        {
            this.expectedValueCutoff = OptionalDouble.of(expectedValueCutoff);
            return this;
        }

        public Builder withWordSize(int wordSize)
        {
            this.wordSize = OptionalInt.of(wordSize);
            return this;
        }

        public Builder withReward(int reward)
        {
            Validate.isTrue(reward >= 0);
            this.reward = OptionalInt.of(reward);
            return this;
        }

        public Builder withPenalty(int penalty)
        {
            Validate.isTrue(penalty <= 0);
            this.penalty = OptionalInt.of(penalty);
            return this;
        }

        public Builder withGapOpen(int gapopen)
        {
            Validate.isTrue(gapopen >= 0);
            this.gapopen = OptionalInt.of(gapopen);
            return this;
        }

        public Builder withGapExtend(int gapextend)
        {
            Validate.isTrue(gapextend >= 0);
            this.gapextend = OptionalInt.of(gapextend);
            return this;
        }

        public Builder withSubjectBestHit(boolean subjectBestHit)
        {
            this.subjectBestHit = subjectBestHit;
            return this;
        }

        public BlastnRunner build()
        {
            if(prefix == null || prefix.isEmpty())
            {
                throw new IllegalStateException("Mandatory field 'prefix' is not set");
            }
            if(blastDir == null || blastDir.isEmpty())
            {
                throw new IllegalStateException("Mandatory field 'blastDir' is not set");
            }
            if(blastDb == null || blastDb.isEmpty())
            {
                throw new IllegalStateException("Mandatory field 'blastDb' is not set");
            }
            if(outputDir == null || outputDir.isEmpty())
            {
                throw new IllegalStateException("Mandatory field 'outputDir' is not set");
            }
            if(numThreads <= 0)
            {
                throw new IllegalStateException("Mandatory field 'numThreads' is not set or invalid");
            }
            return new BlastnRunner(this);
        }
    }
}
