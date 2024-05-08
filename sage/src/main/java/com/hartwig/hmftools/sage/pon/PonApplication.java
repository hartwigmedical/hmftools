package com.hartwig.hmftools.sage.pon;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class PonApplication implements AutoCloseable
{
    private static final String IN_VCF = "in";
    private static final String OUT_VCF = "out";
    private static final String GLOB = "*.sage.somatic.vcf.gz";

    public static void main(String[] args) throws IOException, ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addPath(IN_VCF, true, "Input file");
        configBuilder.addConfigItem(OUT_VCF, true, "Output file");
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        int threads = parseThreads(configBuilder);

        String inputFilePath = configBuilder.getValue(IN_VCF);
        String outputFilePath = configBuilder.getValue(OUT_VCF);

        if(outputFilePath == null || inputFilePath == null)
        {
            SG_LOGGER.error("missing input or output VCFs");
            System.exit(1);
        }

        try(PonApplication app = new PonApplication(threads, inputFilePath, outputFilePath))
        {
            app.run();
        }
    }

    private final PonVCF vcf;
    private final String input;
    private final List<File> files;
    private final ExecutorService executorService;

    private PonApplication(int threads, final String input, final String output) throws IOException
    {
        SG_LOGGER.info("Input: {}", input);
        SG_LOGGER.info("Output: {}", output);

        executorService = Executors.newFixedThreadPool(threads);
        this.input = input;

        files = Lists.newArrayList();
        for(Path path : Files.newDirectoryStream(new File(input).toPath(), GLOB))
        {
            files.add(path.toFile());
        }

        this.vcf = new PonVCF(output, files.size());
    }

    private void run() throws IOException, ExecutionException, InterruptedException
    {
        if(files.isEmpty())
        {
            return;
        }

        final VCFFileReader dictionaryReader = new VCFFileReader(files.get(0), true);
        SAMSequenceDictionary dictionary = dictionaryReader.getFileHeader().getSequenceDictionary();
        dictionaryReader.close();

        /*
        for(SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            SG_LOGGER.info("Processing sequence {}", samSequenceRecord.getSequenceName());
            final PonBuilder ponBuilder = new PonBuilder();

            List<Future<?>> contigFutures = Lists.newArrayList();

            for(Path file : Files.newDirectoryStream(new File(input).toPath(), GLOB))
            {
                Runnable runnable = () -> addVariantsFromFileToBuilder(ponBuilder, samSequenceRecord, file);
                contigFutures.add(executorService.submit(runnableTaskCompletion.task(runnable)));
            }

            for(Future<?> contigFuture : contigFutures)
            {
                contigFuture.get();
            }

            vcf.write(ponBuilder.build());
        }
        */
    }

    private void addVariantsFromFileToBuilder(final PonBuilder ponBuilder, final SAMSequenceRecord samSequenceRecord, final Path file)
    {
        try(VCFFileReader fileReader = new VCFFileReader(file.toFile(), true))
        {
            CloseableIterator<VariantContext> iter =
                    fileReader.query(samSequenceRecord.getSequenceName(), 1, samSequenceRecord.getSequenceLength());
            while(iter.hasNext())
            {
                ponBuilder.add(iter.next());
            }
            iter.close();
        }
    }

    @Override
    public void close()
    {
        executorService.shutdown();
        vcf.close();
        SG_LOGGER.info("PON complete");
    }
}
