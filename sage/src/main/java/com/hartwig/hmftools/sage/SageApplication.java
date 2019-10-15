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
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.ContigContext;

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

    private SageApplication(final Options options, final String... args) throws IOException, ParseException {

        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(cmd);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        vcf = new SageVCF(refGenome, config.outputFile(), config.reference(), config.tumor());

    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        long timeStamp = System.currentTimeMillis();
        final List<ContigContext> contigContexts = Lists.newArrayList();

        SAMSequenceDictionary dictionary = dictionary();
        for (final SAMSequenceRecord samSequenceRecord : dictionary.getSequences()) {
            final String contig = samSequenceRecord.getSequenceName();
            if (HumanChromosome.contains(contig)) {
                int maxPosition = samSequenceRecord.getSequenceLength();
                contigContexts.add(runChromosome(contig, config.regionSliceSize(), maxPosition));
            }
        }

        //        contigContexts.add(runChromosome("17", config.regionSliceSize(), 9000000));
        //        contigContexts.add(runChromosome("17", config.regionSliceSize(), dictionary().getSequence("17").getSequenceLength()));
        //        contigCon`texts.add(runSingleRegion("17", 32000000, 33000000));
        //        contigContexts.add(runSingleRegion("17", 33000001, 34000000));

        for (final ContigContext contigContext : contigContexts) {
            contigContext.write(vcf);
        }

        long timeTaken = System.currentTimeMillis() - timeStamp;
        System.out.println(" in " + timeTaken);
    }

    private SAMSequenceDictionary dictionary() throws IOException {
        SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(config.referenceBam()));
        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

    @NotNull
    private ContigContext runChromosome(@NotNull final String chromosome, int regionSliceSize, int maxPosition) {

        final ContigContext contigContext = new ContigContext(chromosome);
        for (int i = 0; ; i++) {
            int start = 1 + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize, maxPosition);
            contigContext.add(runRegion(chromosome, start, end));

            if (end >= maxPosition) {
                break;
            }
        }

        return contigContext;
    }

    @NotNull
    private ContigContext runSingleRegion(@NotNull final String chromosome, int start, int end) {

        final ContigContext contigContext = new ContigContext(chromosome);
        contigContext.add(runRegion(chromosome, start, end));
        return contigContext;

    }

    @NotNull
    private Future<List<List<AltContext>>> runRegion(@NotNull final String chromosome, int start, int end) {
        final GenomeRegion region = GenomeRegions.create(chromosome, start, end);
        final SagePipeline regionPipeline = new SagePipeline(region, config, executorService, refGenome);
        return regionPipeline.submit();

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
