package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.hotspot.ModifiableVariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;

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

    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        long timeStamp = System.currentTimeMillis();

        // Note: Turns out you need one samreaderfactory per thread!

        // TODO: add referenceSequence(refGenomeFile)
        //        SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(config.tumorBam().get(0)));

        //        GenomeRegion region1 = GenomeRegions.create("17", 1, 13_000_000);
        //        GenomeRegion region2 = GenomeRegions.create("17", 13_000_000, 26_000_000);
        //        GenomeRegion region3 = GenomeRegions.create("17", 26_000_000, 39_000_000);
        //        GenomeRegion region4 = GenomeRegions.create("17", 39_000_000, 52_000_000);
        //        GenomeRegion region5 = GenomeRegions.create("17", 52_000_000, 65_000_000);
        //        GenomeRegion region6 = GenomeRegions.create("17", 65_000_000, 81_195_210);

//        GenomeRegion region1 = GenomeRegions.create("17", 245733, 245733);
        GenomeRegion region1 = GenomeRegions.create("17", 1, 1_000_000);
        GenomeRegion region2 = GenomeRegions.create("17", 1_000_001, 2_000_000);
        GenomeRegion region3 = GenomeRegions.create("17", 2_000_001, 3_000_000);
        GenomeRegion region4 = GenomeRegions.create("17", 3_000_001, 4_000_000);
        GenomeRegion region5 = GenomeRegions.create("17", 4_000_001, 5_000_000);
        GenomeRegion region6 = GenomeRegions.create("17", 5_000_001, 6_000_000);

        SageSamConsumer samConsumer1 = new SageSamConsumer(region1, refGenome);
        SageSamConsumer samConsumer2 = new SageSamConsumer(region2, refGenome);
        SageSamConsumer samConsumer3 = new SageSamConsumer(region3, refGenome);
        SageSamConsumer samConsumer4 = new SageSamConsumer(region4, refGenome);
        SageSamConsumer samConsumer5 = new SageSamConsumer(region5, refGenome);
        SageSamConsumer samConsumer6 = new SageSamConsumer(region6, refGenome);

        List<Future<SageSamConsumer>> futures = Lists.newArrayList();
        //        futures.add(executorService.submit(() -> callable(GenomeRegions.create("5", 1, 60_000_000),  samConsumer)));
        //        futures.add(executorService.submit(() -> callable(GenomeRegions.create("5", 1, 60_000_000),  samConsumer)));
        futures.add(executorService.submit(() -> callable(region1, samConsumer1)));
        futures.add(executorService.submit(() -> callable(region2, samConsumer2)));
        futures.add(executorService.submit(() -> callable(region3, samConsumer3)));
        futures.add(executorService.submit(() -> callable(region4, samConsumer4)));
        futures.add(executorService.submit(() -> callable(region5, samConsumer5)));
        futures.add(executorService.submit(() -> callable(region6, samConsumer6)));
        //        futures.add(executorService.submit(() -> callable(GenomeRegions.create("5", 1, 180915260L), samConsumer)));

//        SageVCF vcf = new SageVCF("/Users/jon/hmf/tmp/colo829.sage.vcf", "NO_NORMAL_YET", "COLO829T");
        for (Future<SageSamConsumer> future : futures) {
            SageSamConsumer consumer = future.get();
//            consumer.evidence().forEach(vcf::write);
        }

//        vcf.close();

        long timeTaken = System.currentTimeMillis() - timeStamp;
        long total = samConsumer1.count() + samConsumer2.count() + samConsumer3.count() + samConsumer4.count() + samConsumer5.count()
                + samConsumer6.count();
        System.out.println(total + " in " + timeTaken);

    }

    private String toString(ModifiableVariantHotspotEvidence evidence) {
        return evidence.chromosome() + "\t" + evidence.position() + "\t" + evidence.ref() + "\t" + evidence.alt() + "\t"
                + evidence.altSupport() + "\t" + evidence.altQuality() + "\t" + evidence.refSupport() + "\t" + evidence.indelSupport();
    }

    private SageSamConsumer callable(GenomeRegion region, SageSamConsumer consumer) throws IOException {
        SamReader tumorReader =
                SamReaderFactory.makeDefault().referenceSequence(new File(config.refGenome())).open(new File(config.tumorBam().get(0)));

        SAMSlicer slicer = new SAMSlicer(13, Lists.newArrayList(region));
        slicer.slice(tumorReader, consumer);
        tumorReader.close();
        return consumer;
    }

    @Override
    public void close() throws IOException {
        refGenome.close();
        executorService.shutdown();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
