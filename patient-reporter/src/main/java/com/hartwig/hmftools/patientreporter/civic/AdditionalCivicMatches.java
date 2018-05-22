package com.hartwig.hmftools.patientreporter.civic;

import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Collection;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AdditionalCivicMatches {
    private static final Logger LOGGER = LogManager.getLogger(AdditionalCivicMatches.class);
    private static final InputStream ADDITIONAL_MATCHES_CSV =
            AdditionalCivicMatches.class.getResourceAsStream("/civic_additional_matches.csv");
    private static final Pattern variantIdPattern = Pattern.compile("/variants/([0-9]+)/");
    private static final String COPY_GAIN = "copy-gain";
    private static final String COPY_LOSS = "copy-loss";
    private static final String LOSS_OF_FUNCTION = "loss-of-function";
    private static final Multimap<String, Integer> ADDITIONAL_VARIANTS_MAPPING;

    static {
        ADDITIONAL_VARIANTS_MAPPING = createAdditionalVariantsMapping();
    }

    @NotNull
    @VisibleForTesting
    static Multimap<String, Integer> createAdditionalVariantsMapping() {
        Multimap<String, Integer> additionalVariantsMapping = ArrayListMultimap.create();
        try {
            final CSVParser parser = CSVParser.parse(ADDITIONAL_MATCHES_CSV, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
            for (final CSVRecord record : parser) {
                final String gene = record.get("gene").toLowerCase();
                final List<Integer> copyGainVariants = variantIdsFromCsvEntry(record.get(COPY_GAIN));
                final List<Integer> copyFullLossVariants = variantIdsFromCsvEntry(record.get(COPY_LOSS));
                final List<Integer> lossOfFunctionVariants = variantIdsFromCsvEntry(record.get(LOSS_OF_FUNCTION));
                additionalVariantsMapping.putAll(gene + COPY_GAIN, copyGainVariants);
                additionalVariantsMapping.putAll(gene + COPY_LOSS, copyFullLossVariants);
                additionalVariantsMapping.putAll(gene + LOSS_OF_FUNCTION, lossOfFunctionVariants);
            }
        } catch (IOException e) {
            LOGGER.error("Could not read civic variant mapping file for copy gain/loss and loss of function. Error message: {}",
                    e.getMessage());
            throw new RuntimeException(e.getMessage());
        }
        return additionalVariantsMapping;
    }

    @NotNull
    private static List<Integer> variantIdsFromCsvEntry(@NotNull final String csvEntry) {
        final List<Integer> variantIds = Lists.newArrayList();
        final Matcher matcher = variantIdPattern.matcher(csvEntry);
        while (matcher.find()) {
            variantIds.add(Integer.parseInt(matcher.group(1)));
        }
        return variantIds;
    }

    private AdditionalCivicMatches() {
    }

    @NotNull
    public static Collection<Integer> copyNumberVariants(@NotNull final GeneCopyNumber copyNumberReport) {
        if (copyNumberReport.alteration() == CopyNumberAlteration.COPY_FULL_LOSS
                || copyNumberReport.alteration() == CopyNumberAlteration.COPY_PARTIAL_LOSS) {
            return copyLossVariants(copyNumberReport.gene());
        } else if (copyNumberReport.alteration() == CopyNumberAlteration.GAIN) {
            return copyGainVariants(copyNumberReport.gene());
        }
        return Lists.newArrayList();
    }

    @NotNull
    public static Collection<Integer> lossOfFunctionVariants(@NotNull final VariantReport variantReport) {
        if (variantReport.consequence().contains(FRAMESHIFT_VARIANT.readableSequenceOntologyTerm()) || variantReport.consequence()
                .contains(STOP_GAINED.readableSequenceOntologyTerm())) {
            return lossOfFunctionVariants(variantReport.gene());
        }
        return Lists.newArrayList();
    }

    @NotNull
    private static Collection<Integer> copyGainVariants(@NotNull final String gene) {
        return ADDITIONAL_VARIANTS_MAPPING.get(gene.toLowerCase() + COPY_GAIN);
    }

    @NotNull
    private static Collection<Integer> copyLossVariants(@NotNull final String gene) {
        return ADDITIONAL_VARIANTS_MAPPING.get(gene.toLowerCase() + COPY_LOSS);
    }

    @NotNull
    private static Collection<Integer> lossOfFunctionVariants(@NotNull final String gene) {
        return ADDITIONAL_VARIANTS_MAPPING.get(gene.toLowerCase() + LOSS_OF_FUNCTION);
    }
}
