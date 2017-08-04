package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.common.variant.VariantFactory.VCF_COLUMN_SEPARATOR;
import static com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariant.ALLELE_SEPARATOR;
import static com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariantFactory.FORMAT_SEPARATOR;
import static com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariantFactory.FORMAT_VALUES_SEPARATOR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariant;
import com.hartwig.hmftools.common.variant.vcf.StrelkaVCFSomaticFile;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class StrelkaPostProcessApplication {

    private static final Logger LOGGER = LogManager.getLogger(StrelkaPostProcessApplication.class);

    private static final String HIGH_CONFIDENCE_BED = "hc_bed";
    private static final String INPUT_VCF = "v";
    private static final String OUTPUT_VCF = "o";
    private static final String SAMPLE_NAME = "t";

    private static final String DP_FIELD_KEY = "DP";
    private static final String GT_VALUE = "0/1";
    private static final String SIMPLIFIED_FORMAT = "GT:AD:DP";
    private static final String AD_FORMAT_HEADER_LINE =
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">";
    private static final String GATK_COMPAT_HEADER_LINE =
            "##StrelkaGATKCompatibility=Added GT fields to strelka calls for gatk compatibility.";
    private static final double THRESHOLD = 1.3;
    private static final int LC_QUALITY_SCORE_THRESHOLD = 20;
    private static final double LC_ALLELE_FREQUENCY_THRESHOLD = 0.1;

    public static void main(final String... args) throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String inputVcf = cmd.getOptionValue(INPUT_VCF);
        final String outputVcf = cmd.getOptionValue(OUTPUT_VCF);
        final String sampleName = cmd.getOptionValue(SAMPLE_NAME);

        if (highConfidenceBed == null || inputVcf == null || outputVcf == null || sampleName == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Strelka Post Process", options);
            System.exit(1);
        }
        final Slicer highConfidenceSlicer = SlicerFactory.fromBedFile(highConfidenceBed);
        processVariants(inputVcf, highConfidenceSlicer, outputVcf, sampleName);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path towards the high confidence bed");
        options.addOption(INPUT_VCF, true, "Path towards the input VCF");
        options.addOption(OUTPUT_VCF, true, "Path towards the output VCF");
        options.addOption(SAMPLE_NAME, true, "Name of the sample");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void processVariants(@NotNull final String filePath, @NotNull final Slicer highConfidenceSlicer,
            @NotNull final String outputVcf, @NotNull final String sampleName) throws IOException, HartwigException {
        final StrelkaVCFSomaticFile strelkaFile = VCFFileLoader.loadStrelkaVCF(filePath);
        final List<String> linesToWrite = Lists.newArrayList();
        final List<String> newVariantLines = strelkaFile.variants()
                .stream()
                .filter(variant -> checkVariant(variant, highConfidenceSlicer))
                .map(StrelkaPostProcessApplication::replaceVariantFormatAndData)
                .collect(Collectors.toList());
        linesToWrite.addAll(strelkaFile.originalMetaInformationLines());
        linesToWrite.add(AD_FORMAT_HEADER_LINE);
        linesToWrite.add(GATK_COMPAT_HEADER_LINE);
        linesToWrite.add(transformHeader(strelkaFile.originalHeaderLine(), sampleName));
        linesToWrite.addAll(newVariantLines);
        Files.write(new File(outputVcf).toPath(), linesToWrite);
        LOGGER.info("Written  output variants to " + outputVcf);
    }

    @VisibleForTesting
    static boolean checkVariant(@NotNull final StrelkaSomaticVariant variant, @NotNull final Slicer highConfidenceSlicer) {
        if (variant.alt().split(ALLELE_SEPARATOR).length > 1) {
            LOGGER.warn("More than 1 alt for record: " + variant.originalVCFLine());
            return true;
        } else if (variant.alt().equals(".")) {
            LOGGER.warn("Alt is . for record: " + variant.originalVCFLine());
            return false;
        }
        try {
            return variant.qualityScore() > LC_QUALITY_SCORE_THRESHOLD && variant.allelicFrequency() > LC_ALLELE_FREQUENCY_THRESHOLD || (
                    highConfidenceSlicer.includes(variant) && variant.allelicFrequency() * variant.qualityScore() > THRESHOLD);
        } catch (final HartwigException e) {
            LOGGER.error(e.getMessage());
            return false;
        }
    }

    @NotNull
    @VisibleForTesting
    static String simplifyTumorDataString(@NotNull final StrelkaSomaticVariant variant) {
        final String dp = variant.tumorData().get(DP_FIELD_KEY).get(0);
        final String ad_ref = variant.readRefAD();
        final String ad_alt = variant.readAltAD();
        return GT_VALUE + FORMAT_SEPARATOR + ad_ref + FORMAT_VALUES_SEPARATOR + ad_alt + FORMAT_SEPARATOR + dp;
    }

    @NotNull
    private static String replaceVariantFormatAndData(@NotNull final StrelkaSomaticVariant variant) {
        final List<String> vcfFields = Lists.newArrayList(variant.originalVCFLine().split(VCF_COLUMN_SEPARATOR));
        // MIVO: remove format, normal & tumor fields
        vcfFields.remove(vcfFields.size() - 1);
        vcfFields.remove(vcfFields.size() - 1);
        vcfFields.remove(vcfFields.size() - 1);
        vcfFields.add(SIMPLIFIED_FORMAT);
        vcfFields.add(simplifyTumorDataString(variant));
        return Strings.join(vcfFields, VCF_COLUMN_SEPARATOR.charAt(0));
    }

    @NotNull
    private static String transformHeader(@NotNull final String headerLine, @NotNull final String sampleName) {
        final List<String> headerFields = Lists.newArrayList(headerLine.split(VCF_COLUMN_SEPARATOR));
        // MIVO: remove normal & tumor fields
        headerFields.remove(headerFields.size() - 1);
        headerFields.remove(headerFields.size() - 1);
        headerFields.add(sampleName);
        return Strings.join(headerFields, VCF_COLUMN_SEPARATOR.charAt(0));
    }
}
