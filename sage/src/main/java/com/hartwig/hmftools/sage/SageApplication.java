package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.task.SagePipeline;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final SageVCF vcf;

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException {
        final Options options = SageConfig.createOptions();
        try (final SageApplication application = new SageApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("AmberApplication", options);
            System.exit(1);
        }
    }

    public SageApplication(final Options options, final String... args) throws IOException, ParseException {

        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(cmd);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        vcf = new SageVCF(refGenome, "/Users/jon/hmf/tmp/colo829.sage.vcf", config.reference(), config.tumor());

    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        long timeStamp = System.currentTimeMillis();

        // Note: Turns out you need one samreaderfactory per thread!

        List<Future<List<List<AltContext>>>> futures = Lists.newArrayList();

        for (int j = 0; j < 6; j++) {
            int start = 1 + j * 1_000_000;
            int end = 1_000_000 + j * 1_000_000;

            final GenomeRegion region = GenomeRegions.create("17", start, end);
            SagePipeline myThing = new SagePipeline(region, config, executorService, refGenome);

            futures.add(myThing.submit());
        }

        for (Future<List<List<AltContext>>> future : futures) {
            //            future.get().forEach(System.out::println);
            future.get().forEach(vcf::write);
        }

        long timeTaken = System.currentTimeMillis() - timeStamp;
        System.out.println(" in " + timeTaken);
    }

    private <T extends Consumer<SAMRecord>> T callable(GenomeRegion region, T consumer, String bamFile) throws IOException {
        SamReader tumorReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.refGenome())).open(new File(bamFile));

        SAMSlicer slicer = new SAMSlicer(0, Lists.newArrayList(region));
        slicer.slice(tumorReader, consumer);
        tumorReader.close();
        return consumer;
    }

    private <T extends Consumer<SAMRecord>> T callable(List<GenomeRegion> regions, T consumer, String bamFile) throws IOException {
        SamReader tumorReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.refGenome())).open(new File(bamFile));

        SAMSlicer slicer = new SAMSlicer(13, regions);
        slicer.slice(tumorReader, consumer);
        tumorReader.close();
        return consumer;
    }

    @Override
    public void close() throws IOException {
        vcf.close();
        refGenome.close();
        executorService.shutdown();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
