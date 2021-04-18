package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.sage.SageApplication.createCommandLine;

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
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.pipeline.AdditionalReferencePipeline;
import com.hartwig.hmftools.sage.pipeline.ChromosomePartition;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationSupplier;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.apache.commons.cli.CommandLine;
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
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class SageAppendApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageAppendApplication.class);
    private static final double MIN_PRIOR_VERSION = 2.4;

    public static void main(String[] args) {
        final Options options = SageConfig.createAddReferenceOptions();
        try (final SageAppendApplication application = new SageAppendApplication(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageAppendApplication", options);
            System.exit(1);
        } catch (Exception e) {
            LOGGER.warn(e);
            System.exit(1);
        }
    }

    private final SageVCF outputVCF;
    private final SageConfig config;
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final QualityRecalibrationSupplier qualityRecalibrationSupplier;
    private final AbstractFeatureReader<VariantContext, LineIterator> inputReader;
    private final long timeStamp = System.currentTimeMillis();

    public SageAppendApplication(final Options options, final String... args) throws ParseException, IOException {
        final VersionInfo version = new VersionInfo("sage.version");
        LOGGER.info("SAGE version: {}", version.version());

        final CommandLine cmd = createCommandLine(args, options);
        this.config = SageConfig.createConfig(true, version.version(), cmd);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        executorService = Executors.newFixedThreadPool(config.threads(), namedThreadFactory);
        refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        qualityRecalibrationSupplier = new QualityRecalibrationSupplier(executorService, refGenome, config);

        final String inputVcf = config.inputFile();
        inputReader = AbstractFeatureReader.getFeatureReader(inputVcf, new VCFCodec(), false);

        VCFHeader inputHeader = (VCFHeader) inputReader.getHeader();
        LOGGER.info("Reading and validating file: {}", inputVcf);
        validateInputHeader(inputHeader);

        outputVCF = new SageVCF(refGenome, config, inputHeader);
        LOGGER.info("Writing to file: {}", config.outputFile());
    }

    public void run() throws IOException, ExecutionException, InterruptedException {
        final ChromosomePartition chromosomePartition = new ChromosomePartition(config, refGenome);
        final List<VariantContext> existing = verifyAndReadExisting();

        final SAMSequenceDictionary dictionary = dictionary();
        final List<Future<List<VariantContext>>> futures = Lists.newArrayList();
        final Map<String, QualityRecalibrationMap> recalibrationMap = qualityRecalibrationSupplier.get();
        final AdditionalReferencePipeline pipeline = new AdditionalReferencePipeline(config, executorService, refGenome, recalibrationMap);

        for (final SAMSequenceRecord samSequenceRecord : dictionary.getSequences()) {
            final String contig = samSequenceRecord.getSequenceName();
            if (HumanChromosome.contains(contig) || MitochondrialChromosome.contains(contig)) {
                final List<VariantContext> chromosomeVariants =
                        existing.stream().filter(x -> x.getContig().equals(contig)).collect(Collectors.toList());

                for (GenomeRegion region : chromosomePartition.partition(contig)) {
                    final List<VariantContext> regionVariants = chromosomeVariants.stream()
                            .filter(x -> x.getStart() >= region.start() && x.getStart() <= region.end())
                            .collect(Collectors.toList());

                    futures.add(pipeline.appendReference(region, regionVariants));
                }
            }
        }

        for (Future<List<VariantContext>> updatedVariantsFuture : futures) {
            final List<VariantContext> updatedVariants = updatedVariantsFuture.get();
            updatedVariants.forEach(outputVCF::write);
        }
    }

    @Override
    public void close() throws IOException {
        inputReader.close();
        outputVCF.close();
        refGenome.close();
        executorService.shutdown();
        long timeTaken = System.currentTimeMillis() - timeStamp;
        LOGGER.info("Completed in {} seconds", timeTaken / 1000);
    }

    @NotNull
    public List<VariantContext> verifyAndReadExisting() throws IOException, IllegalArgumentException {
        VCFHeader header = (VCFHeader) inputReader.getHeader();

        List<VariantContext> result = Lists.newArrayList();
        for (VariantContext variantContext : inputReader.iterator()) {
            result.add(variantContext.fullyDecode(header, false));
        }

        return result;
    }

    public void validateInputHeader(VCFHeader header) throws IllegalArgumentException {
        double oldVersion = sageVersion(header);
        if (Doubles.lessThan(oldVersion, MIN_PRIOR_VERSION)) {
            throw new IllegalArgumentException("Input VCF must be from SAGE version " + MIN_PRIOR_VERSION + " onwards");
        }

        final Set<String> samplesInExistingVcf = existingSamples(header);
        for (String refSample : config.reference()) {
            if (samplesInExistingVcf.contains(refSample)) {
                throw new IllegalArgumentException("Sample " + refSample + " already exits in input VCF");
            }
        }
    }

    private static double sageVersion(@NotNull final VCFHeader header) {
        VCFHeaderLine oldVersion = header.getMetaDataLine(SageVCF.VERSION_META_DATA);
        if (oldVersion == null) {
            return 0;
        }

        String oldVersionString = oldVersion.getValue();
        try {
            return Double.parseDouble(oldVersionString);
        } catch (Exception e) {
            return 0;
        }
    }

    @NotNull
    private static Set<String> existingSamples(@NotNull final VCFHeader header) {
        return Sets.newHashSet(header.getGenotypeSamples());
    }

    private SAMSequenceDictionary dictionary() throws IOException {
        final String bam = config.referenceBam().isEmpty() ? config.tumorBam().get(0) : config.referenceBam().get(0);
        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(config.validationStringency())
                .referenceSource(new ReferenceSource(refGenome)).open(new File(bam));
        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }

}
