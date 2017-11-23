package com.hartwig.hmftools.patientreporter.civic;

import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReportType.GAIN;
import static com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReportType.LOSS;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Collection;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
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
    private static final Multimap<String, Integer> additionalVariantsMapping;

    static {
        additionalVariantsMapping = ArrayListMultimap.create();
        try {
            final CSVParser parser = CSVParser.parse(ADDITIONAL_MATCHES_CSV, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
            for (final CSVRecord record : parser) {
                final String gene = record.get("gene").toLowerCase();
                final List<Integer> copyGainVariants = variantIdsFromCsvEntry(record.get(COPY_GAIN));
                final List<Integer> copyLossVariants = variantIdsFromCsvEntry(record.get(COPY_LOSS));
                final List<Integer> lossOfFunctionVariants = variantIdsFromCsvEntry(record.get(LOSS_OF_FUNCTION));
                additionalVariantsMapping.putAll(gene + COPY_GAIN, copyGainVariants);
                additionalVariantsMapping.putAll(gene + COPY_LOSS, copyLossVariants);
                additionalVariantsMapping.putAll(gene + LOSS_OF_FUNCTION, lossOfFunctionVariants);
            }
        } catch (IOException e) {
            LOGGER.error("Could not read civic variant mapping file for copy gain/loss and loss of function. Error message: {}",
                    e.getMessage());
            throw new RuntimeException(e.getMessage());
        }
    }

    private AdditionalCivicMatches() {
    }

    private static List<Integer> variantIdsFromCsvEntry(@NotNull final String csvEntry) {
        final List<Integer> variantIds = Lists.newArrayList();
        final Matcher matcher = variantIdPattern.matcher(csvEntry);
        while (matcher.find()) {
            variantIds.add(Integer.parseInt(matcher.group(1)));
        }
        return variantIds;
    }

    public static Collection<Integer> copyNumberVariants(@NotNull final CopyNumberReport copyNumberReport) {
        if (copyNumberReport.type() == LOSS) {
            return copyLossVariants(copyNumberReport.gene());
        } else if (copyNumberReport.type() == GAIN) {
            return copyGainVariants(copyNumberReport.gene());
        }
        return Lists.newArrayList();
    }

    public static Collection<Integer> lossOfFunctionVariants(@NotNull final VariantReport variantReport) {
        if (variantReport.consequence().contains(FRAMESHIFT_VARIANT.readableSequenceOntologyTerm()) || variantReport.consequence()
                .contains(STOP_GAINED.readableSequenceOntologyTerm())) {
            return lossOfFunctionVariants(variantReport.gene());
        }
        return Lists.newArrayList();
    }

    private static Collection<Integer> copyGainVariants(@NotNull final String gene) {
        return additionalVariantsMapping.get(gene.toLowerCase() + COPY_GAIN);
    }

    private static Collection<Integer> copyLossVariants(@NotNull final String gene) {
        return additionalVariantsMapping.get(gene.toLowerCase() + COPY_LOSS);
    }

    private static Collection<Integer> lossOfFunctionVariants(@NotNull final String gene) {
        return additionalVariantsMapping.get(gene.toLowerCase() + LOSS_OF_FUNCTION);
    }
}
