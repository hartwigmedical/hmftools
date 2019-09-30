package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.sage.count.BaseDetails;
import com.hartwig.hmftools.sage.count.ReadContextConsumer;
import com.hartwig.hmftools.sage.count.ReadContextConsumerDispatcher;
import com.hartwig.hmftools.sage.count.ReadContextCount;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
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

        LOGGER.info("Examining tumor sample for evidence of variants");
        final List<BaseDetails> tumorBaseDetails = tumor(config.tumorBam().get(0));

        repeatContextStuff(config.tumorBam().get(0), tumorBaseDetails);

        final Set<Long> hotspots = tumorBaseDetails.stream().map(BaseDetails::position).collect(Collectors.toSet());

        LOGGER.info("Examining normal sample for evidence of variants");
        final List<BaseDetails> normalEvidence = normal(config.referenceBam(), hotspots);

        LOGGER.info("Combining");
        SageVCF vcf = new SageVCF(refGenome, "/Users/jon/hmf/tmp/colo829.sage.vcf", "COLO829R", "COLO829T");

        Map<Long, BaseDetails> normalMap = asMap(normalEvidence);
        for (final BaseDetails tumorBase : tumorBaseDetails) {

            @Nullable
            final BaseDetails normalBase = normalMap.get(tumorBase.position());
            for (VariantHotspotEvidence tumorHotspot : tumorBase.evidence()) {
                @Nullable
                final VariantHotspotEvidence normalHotspot =
                        normalBase == null ? null : normalBase.selectOrCreate(tumorHotspot.ref(), tumorHotspot.alt());
                vcf.write(tumorHotspot, normalHotspot);

            }
        }

        vcf.close();

        long timeTaken = System.currentTimeMillis() - timeStamp;

        System.out.println(" in " + timeTaken);
    }

    @NotNull
    private List<BaseDetails> tumor(String bamFile) throws ExecutionException, InterruptedException {

        GenomeRegion region1 = GenomeRegions.create("17", 1, 1_000_000);
        GenomeRegion region2 = GenomeRegions.create("17", 1_000_001, 2_000_000);
        GenomeRegion region3 = GenomeRegions.create("17", 2_000_001, 3_000_000);
        GenomeRegion region4 = GenomeRegions.create("17", 3_000_001, 4_000_000);
        GenomeRegion region5 = GenomeRegions.create("17", 4_000_001, 5_000_000);
        GenomeRegion region6 = GenomeRegions.create("17", 5_000_001, 6_000_000);

        SageSamConsumer samConsumer1 = new SageSamConsumer(13, region1, refGenome);
        SageSamConsumer samConsumer2 = new SageSamConsumer(13, region2, refGenome);
        SageSamConsumer samConsumer3 = new SageSamConsumer(13, region3, refGenome);
        SageSamConsumer samConsumer4 = new SageSamConsumer(13, region4, refGenome);
        SageSamConsumer samConsumer5 = new SageSamConsumer(13, region5, refGenome);
        SageSamConsumer samConsumer6 = new SageSamConsumer(13, region6, refGenome);

        List<Future<SageSamConsumer>> futures = Lists.newArrayList();
        futures.add(executorService.submit(() -> callable(region1, samConsumer1, bamFile)));
        futures.add(executorService.submit(() -> callable(region2, samConsumer2, bamFile)));
        futures.add(executorService.submit(() -> callable(region3, samConsumer3, bamFile)));
        futures.add(executorService.submit(() -> callable(region4, samConsumer4, bamFile)));
        futures.add(executorService.submit(() -> callable(region5, samConsumer5, bamFile)));
        futures.add(executorService.submit(() -> callable(region6, samConsumer6, bamFile)));

        final List<BaseDetails> tumorEvidence = Lists.newArrayList();

        for (Future<SageSamConsumer> future : futures) {
            SageSamConsumer consumer = future.get();
            tumorEvidence.addAll(consumer.bases());
        }

        return tumorEvidence;
    }

    @NotNull
    private List<BaseDetails> normal(String bamFile, Set<Long> hotspots) throws ExecutionException, InterruptedException {

        GenomeRegion region1 = GenomeRegions.create("17", 1, 1_000_000);
        GenomeRegion region2 = GenomeRegions.create("17", 1_000_001, 2_000_000);
        GenomeRegion region3 = GenomeRegions.create("17", 2_000_001, 3_000_000);
        GenomeRegion region4 = GenomeRegions.create("17", 3_000_001, 4_000_000);
        GenomeRegion region5 = GenomeRegions.create("17", 4_000_001, 5_000_000);
        GenomeRegion region6 = GenomeRegions.create("17", 5_000_001, 6_000_000);

        SageSamConsumer samConsumer1 = new SageSamConsumer(13, region1, refGenome, hotspots);
        SageSamConsumer samConsumer2 = new SageSamConsumer(13, region2, refGenome, hotspots);
        SageSamConsumer samConsumer3 = new SageSamConsumer(13, region3, refGenome, hotspots);
        SageSamConsumer samConsumer4 = new SageSamConsumer(13, region4, refGenome, hotspots);
        SageSamConsumer samConsumer5 = new SageSamConsumer(13, region5, refGenome, hotspots);
        SageSamConsumer samConsumer6 = new SageSamConsumer(13, region6, refGenome, hotspots);

        List<Future<SageSamConsumer>> futures = Lists.newArrayList();
        futures.add(executorService.submit(() -> callable(region1, samConsumer1, bamFile)));
        futures.add(executorService.submit(() -> callable(region2, samConsumer2, bamFile)));
        futures.add(executorService.submit(() -> callable(region3, samConsumer3, bamFile)));
        futures.add(executorService.submit(() -> callable(region4, samConsumer4, bamFile)));
        futures.add(executorService.submit(() -> callable(region5, samConsumer5, bamFile)));
        futures.add(executorService.submit(() -> callable(region6, samConsumer6, bamFile)));

        final List<BaseDetails> tumorEvidence = Lists.newArrayList();

        for (Future<SageSamConsumer> future : futures) {
            SageSamConsumer consumer = future.get();
            tumorEvidence.addAll(consumer.bases());
        }

        return tumorEvidence;
    }

    private void repeatContextStuff(@NotNull String bamFile, List<BaseDetails> tumorDetails)
            throws ExecutionException, InterruptedException {

        long time = System.currentTimeMillis();
        LOGGER.info("Getting repeat contexts...");

        // TODO: This was a terrible idea and took 15 mins....

        List<Future<ReadContextConsumerDispatcher>> futures = Lists.newArrayList();

        for (int j = 0; j < 6; j++) {
            int start = 1 + j * 1_000_000;
            int end = 1_000_000 + j * 1_000_000;
            GenomeRegion region = GenomeRegions.create("17", start, end);

            GenomeRegions regionsBuilder = new GenomeRegions("17", 1000);

            List<ReadContextConsumer> readContexts = Lists.newArrayList();
            for (BaseDetails detail : tumorDetails) {
                for (VariantHotspotEvidence evidence : detail.evidence()) {
                    if (evidence.altSupport() > 2 && region.contains(evidence)) {

                        List<ReadContextCount> altReadContexts = detail.contexts(evidence.alt());
                        if (!altReadContexts.isEmpty()) {
                            regionsBuilder.addPosition(evidence.position());

                            readContexts.add(new ReadContextConsumer(evidence, altReadContexts.get(0)));
                        }
                    }
                }
            }

            ReadContextConsumerDispatcher dispatcher = new ReadContextConsumerDispatcher(readContexts);
            futures.add(executorService.submit(() -> callable(regionsBuilder.build(), dispatcher, bamFile)));
        }

        LOGGER.info("submitted, just awaiting results");
        final List<ReadContextConsumer> result = Lists.newArrayList();
        for (Future<ReadContextConsumerDispatcher> future : futures) {
            result.addAll(future.get().consumers());
        }


        LOGGER.info("Getting repeat contexts complete in {} millis!", System.currentTimeMillis() - time);
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
        refGenome.close();
        executorService.shutdown();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Map<Long, BaseDetails> asMap(@NotNull final List<? extends BaseDetails> evidence) {
        return evidence.stream().collect(Collectors.toMap(BaseDetails::position, x -> x));
    }

}
