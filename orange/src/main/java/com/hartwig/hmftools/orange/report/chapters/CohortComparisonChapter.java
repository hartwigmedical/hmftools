package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.cuppa.CuppaEntry;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CohortComparisonChapter implements ReportChapter {

    private static final Logger LOGGER = LogManager.getLogger(CohortComparisonChapter.class);
    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    private static final String COMBINED_CATEGORY = "COMBINED";
    private static final String CLASSIFIER_CATEGORY = "CLASSIFIER";
    private static final String SAMPLE_TRAIT_CATEGORY = "SAMPLE_TRAIT";
    private static final String FEATURE_CATEGORY = "FEATURE";
    private static final String SV_CATEGORY = "SV";

    @NotNull
    private final OrangeReport report;

    public CohortComparisonChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Cohort Comparison";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Cohort Comparison").addStyle(ReportResources.chapterTitleStyle()));
        Set<String> refCancerTypes = determineRefCancerTypes();

        addTumorTypeClassificationTable(document, refCancerTypes);
        addTumorTraitsTable(document, refCancerTypes);
        addTumorDriverTable(document, refCancerTypes);

        document.add(new Paragraph("TODO: Add detailed cohort incidence per signature").addStyle(ReportResources.tableTitleStyle()));
    }

    private void addTumorTypeClassificationTable(@NotNull Document document, @NotNull Set<String> refCancerTypes) {
        Table table = TableUtil.createReportContentTable(new float[] { 3, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Ref Cancer Type"), TableUtil.createHeaderCell("Combined"),
                        TableUtil.createHeaderCell("SNV 96"), TableUtil.createHeaderCell("Position"), TableUtil.createHeaderCell("Feature"),
                        TableUtil.createHeaderCell("Expression"), TableUtil.createHeaderCell("Alt SJ") });

        for (String refCancerType : refCancerTypes) {
            table.addCell(TableUtil.createContentCell(refCancerType));
            table.addCell(TableUtil.createContentCell(findCombinedLikelihood("DNA_COMBINED", refCancerType)));
            table.addCell(TableUtil.createContentCell(findClassifierLikelihood("SNV_96_PAIRWISE_SIMILARITY", refCancerType)));
            table.addCell(TableUtil.createContentCell(findClassifierLikelihood("GENOMIC_POSITION_SIMILARITY", refCancerType)));
            table.addCell(TableUtil.createContentCell(findClassifierLikelihood("FEATURE", refCancerType)));
            table.addCell(TableUtil.createContentCell("NA"));
            table.addCell(TableUtil.createContentCell("NA"));
        }

        document.add(TableUtil.createWrappingReportTable(table, "Primary tumor classification"));
    }

    private void addTumorTraitsTable(@NotNull Document document, @NotNull Set<String> refCancerTypes) {
        String firstCancerType = Lists.newArrayList(refCancerTypes.iterator()).get(0);
        Table traits = new Table(UnitValue.createPercentArray(new float[] { 2, 1, 3 }));
        traits.addCell(TableUtil.createKeyCell("Gender"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SAMPLE_TRAIT_CATEGORY,
                Sets.newHashSet("GENDER"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("SNVs"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SAMPLE_TRAIT_CATEGORY,
                Sets.newHashSet("SNV_COUNT"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("MSI"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SAMPLE_TRAIT_CATEGORY,
                Sets.newHashSet("MS_INDELS_TMB"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("DUPs 32B-200B"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SV_CATEGORY,
                Sets.newHashSet("SIMPLE_DUP_32B_200B"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("Max Complex Size"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SV_CATEGORY,
                Sets.newHashSet("MAX_COMPLEX_SIZE"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("LINEs"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SV_CATEGORY, Sets.newHashSet("LINE"), null, firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        traits.addCell(TableUtil.createKeyCell("Telomeric SGLs"));
        traits.addCell(TableUtil.createValueCell(safeValue(findEntry(SV_CATEGORY,
                Sets.newHashSet("TELOMERIC_SGL"),
                null,
                firstCancerType))));
        traits.addCell(TableUtil.createValueCell(""));
        document.add(TableUtil.createWrappingReportTable(traits, "Sample traits"));

        Table percentiles = TableUtil.createReportContentTable(new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Ref Cancer Type"), TableUtil.createHeaderCell("Gender"),
                        TableUtil.createHeaderCell("SNV"), TableUtil.createHeaderCell("MSI"), TableUtil.createHeaderCell("DUP"),
                        TableUtil.createHeaderCell("Complex"), TableUtil.createHeaderCell("LINE"),
                        TableUtil.createHeaderCell("Telomeric") });

        for (String refCancerType : refCancerTypes) {
            percentiles.addCell(TableUtil.createContentCell(refCancerType));
            percentiles.addCell(TableUtil.createContentCell(findSampleTraitValue("GENDER", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSampleTraitValue("SNV_COUNT", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSampleTraitValue("MS_INDELS_TMB", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSvValue("SIMPLE_DUP_32B_200B", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSvValue("MAX_COMPLEX_SIZE", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSvValue("LINE", refCancerType)));
            percentiles.addCell(TableUtil.createContentCell(findSvValue("TELOMERIC_SGL", refCancerType)));
        }

        document.add(TableUtil.createWrappingReportTable(percentiles, "Sample trait percentiles"));
    }

    private void addTumorDriverTable(@NotNull Document document, @NotNull Set<String> refCancerTypes) {
        Set<String> drivers = Sets.newHashSet();
        for (CuppaEntry entry : report.cuppaEntries()) {
            if (entry.category().equals(FEATURE_CATEGORY) && (entry.dataType().equals("DRIVER") || entry.dataType().equals("INDEL"))) {
                drivers.add(entry.value());
            }
        }

        List<String> driverList = Lists.newArrayList(drivers.iterator());
        int driversPerTable = 5;
        for (int i = 0; i <= driverList.size() / driversPerTable; i++) {
            renderDriverPrevalenceTable(document,
                    refCancerTypes,
                    driversPerTable,
                    driverList.subList(i * driversPerTable, Math.min((i + 1) * driversPerTable, driverList.size())));
        }
    }

    private void renderDriverPrevalenceTable(@NotNull Document document, @NotNull Set<String> refCancerTypes, int requiredDriverColumns,
            @NotNull List<String> drivers) {
        float[] widths = new float[requiredDriverColumns + 1];
        Cell[] headers = new Cell[requiredDriverColumns + 1];

        widths[0] = 3;
        headers[0] = TableUtil.createHeaderCell("Ref Cancer Type");
        for (int i = 0; i < requiredDriverColumns; i++) {
            widths[i + 1] = 1;
            headers[i + 1] = TableUtil.createHeaderCell(i < drivers.size() ? drivers.get(i) : Strings.EMPTY);
        }

        Table percentiles = TableUtil.createReportContentTable(widths, headers);

        for (String refCancerType : refCancerTypes) {
            percentiles.addCell(TableUtil.createContentCell(refCancerType));
            for (int i = 0; i < requiredDriverColumns; i++) {
                if (i < drivers.size()) {
                    percentiles.addCell(TableUtil.createContentCell(findDriverFeature(drivers.get(i), refCancerType)));
                } else {
                    percentiles.addCell(TableUtil.createContentCell(Strings.EMPTY));
                }
            }
        }

        document.add(TableUtil.createWrappingReportTable(percentiles, "Sample driver prevalences"));
    }

    @NotNull
    private String findDriverFeature(@NotNull String driver, @NotNull String refCancerType) {
        return findPercentileOrPrevalence(FEATURE_CATEGORY, Sets.newHashSet("DRIVER", "INDEL"), driver, refCancerType);
    }

    @NotNull
    private String findCombinedLikelihood(@NotNull String type, @NotNull String refCancerType) {
        return findLikelihood(COMBINED_CATEGORY, type, refCancerType);
    }

    @NotNull
    private String findClassifierLikelihood(@NotNull String type, @NotNull String refCancerType) {
        return findLikelihood(CLASSIFIER_CATEGORY, type, refCancerType);
    }

    @NotNull
    private String findSampleTraitValue(@NotNull String type, @NotNull String refCancerType) {
        return findPercentileOrPrevalence(SAMPLE_TRAIT_CATEGORY, Sets.newHashSet(type), null, refCancerType);
    }

    @NotNull
    private String findSvValue(@NotNull String type, @NotNull String refCancerType) {
        return findPercentileOrPrevalence(SV_CATEGORY, Sets.newHashSet(type), null, refCancerType);
    }

    @NotNull
    private String findLikelihood(@NotNull String category, @NotNull String type, @NotNull String refCancerType) {
        CuppaEntry entry = findEntry(category, Sets.newHashSet(type), null, refCancerType);
        if (entry == null) {
            return "NA";
        } else {
            if (!entry.resultType().equals("LIKELIHOOD")) {
                LOGGER.warn("Not a likelihood entry: {}", entry);
            }
            return PERCENTAGE_FORMAT.format(entry.refValue() * 100);
        }
    }

    @NotNull
    private String findPercentileOrPrevalence(@NotNull String category, @NotNull Set<String> types, @Nullable String value,
            @NotNull String refCancerType) {
        CuppaEntry entry = findEntry(category, types, value, refCancerType);
        if (entry == null) {
            return "NA";
        } else {
            if (!entry.resultType().equals("PERCENTILE") && !entry.resultType().equals("PREVALENCE")) {
                LOGGER.warn("Not a percentile or prevalence entry: {}", entry);
            }
            return PERCENTAGE_FORMAT.format(entry.refValue() * 100);
        }
    }

    @Nullable
    private CuppaEntry findEntry(@NotNull String category, @NotNull Set<String> dataTypes, @Nullable String value,
            @NotNull String refCancerType) {
        for (CuppaEntry entry : report.cuppaEntries()) {
            if (entry.category().equals(category) && dataTypes.contains(entry.dataType()) && entry.refCancerType().equals(refCancerType)
                    && (value == null || entry.value().equals(value))) {
                return entry;
            }
        }
        LOGGER.warn("Could not find entry of category '{}' and datatypes '{}' with value '{}' for ref cancer type '{}'",
                category,
                dataTypes,
                value,
                refCancerType);

        return null;
    }

    @NotNull
    private static String safeValue(@NotNull CuppaEntry entry) {
        return entry != null ? entry.value() : "NA";
    }

    @NotNull
    private Set<String> determineRefCancerTypes() {
        Set<String> refCancerTypes = Sets.newTreeSet(Comparator.naturalOrder());
        for (CuppaEntry data : report.cuppaEntries()) {
            refCancerTypes.add(data.refCancerType());
        }
        return refCancerTypes;
    }
}
