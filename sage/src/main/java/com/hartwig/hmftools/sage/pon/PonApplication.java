package com.hartwig.hmftools.sage.pon;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

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

    public static void main(String[] args) throws IOException, ParseException, ExecutionException, InterruptedException
    {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(IN_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);
        final int threads = parseThreads(cmd, 5);

        if(outputFilePath == null || inputFilePath == null)
        {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PonApplication", options);
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

    private PonApplication(int threads, @NotNull final String input, @NotNull final String output) throws IOException
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

        for(SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            SG_LOGGER.info("Processing sequence {}", samSequenceRecord.getSequenceName());
            final PonBuilder ponBuilder = new PonBuilder();
            final RunnableTaskCompletion runnableTaskCompletion = new RunnableTaskCompletion();

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

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(IN_VCF, true, "Input file");
        options.addOption(OUT_VCF, true, "Output file");
        addThreadOptions(options);
        return options;
    }

    @Override
    public void close()
    {
        executorService.shutdown();
        vcf.close();
        SG_LOGGER.info("PON complete");
    }
}
