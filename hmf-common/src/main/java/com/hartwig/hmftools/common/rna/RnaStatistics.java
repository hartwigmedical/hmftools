package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.FAIL_LOW_COVERAGE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.PASS;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.WARN_DUPLICATE_RATE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.WARN_SPLICED_GENE_COVERAGE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersFromString;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersToString;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferHeaderDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getDoubleValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getLongValue;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.immutables.value.Value;

@Value.Immutable
public abstract class RnaStatistics
{
    public abstract long totalFragments();
    public abstract long duplicateFragments();
    public abstract double splicedFragmentPerc();
    public abstract double unsplicedFragmentPerc();
    public abstract double altFragmentPerc();
    public abstract double chimericFragmentPerc();
    public abstract int splicedGeneCount();

    public abstract int readLength();

    // 5th, 50th and 95th percentile intronic fragment length
    public abstract double fragmentLength5thPercent();
    public abstract double fragmentLength50thPercent();
    public abstract double fragmentLength95thPercent();

    // proportion of fragments in 7 highly expressed genes
    public abstract double enrichedGenePercent();

    // Median GC (excluding 7 highly expressed genes)
    public abstract double medianGCRatio();

    public abstract double forwardStrandPercent();

    public abstract List<RnaQcFilter> qcStatus();

    public static final int LOW_COVERAGE_THRESHOLD = 2500000;
    public static final int SPLICE_GENE_THRESHOLD = 17000;
    private static final double HIGH_DUPLICATES_THRESHOLD = 0.9;

    public static final String SUMMARY_FILE_ID = "summary.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + SUMMARY_FILE_ID;
    }

    public static String csvHeader()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("QcStatus")
                .add("TotalFragments")
                .add("DuplicateFragments")
                .add("SplicedFragmentPerc")
                .add("UnsplicedFragmentPerc")
                .add("AltFragmentPerc")
                .add("ChimericFragmentPerc")
                .add("SplicedGeneCount")
                .add("ReadLength")
                .add("FragLength5th")
                .add("FragLength50th")
                .add("FragLength95th")
                .add("EnrichedGenePercent")
                .add("MedianGCRatio")
                .add("ForwardStrandPercent")
                .toString();
    }

    public String toCsv(final String sampleId)
    {
        return new StringJoiner(DELIMITER)
                .add(sampleId)
                .add(qcFiltersToString(qcStatus()))
                .add(String.valueOf(totalFragments()))
                .add(String.valueOf(duplicateFragments()))
                .add(String.format("%.3f", splicedFragmentPerc()))
                .add(String.format("%.3f", unsplicedFragmentPerc()))
                .add(String.format("%.3f", altFragmentPerc()))
                .add(String.format("%.3f", chimericFragmentPerc()))
                .add(String.valueOf(splicedGeneCount()))
                .add(String.valueOf(readLength()))
                .add(String.format("%.0f", fragmentLength5thPercent()))
                .add(String.format("%.0f", fragmentLength50thPercent()))
                .add(String.format("%.0f", fragmentLength95thPercent()))
                .add(String.format("%.3f", enrichedGenePercent()))
                .add(String.format("%.3f", medianGCRatio()))
                .add(String.format("%.3f", forwardStrandPercent()))
                .toString();
    }

    public static RnaStatistics fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        String fileDelim = inferHeaderDelimiter(header);
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

        final String[] values = lines.get(1).split(fileDelim);

        long totalFragments = getLongValue(fieldsIndexMap, "TotalFragments", values);
        long duplicateFragments = getLongValue(fieldsIndexMap, "DuplicateFragments", values);

        int splicedGeneCount = fieldsIndexMap.containsKey("SplicedGeneCount") ?
                getIntValue(fieldsIndexMap, "SplicedGeneCount", values) : -1;

        List<RnaQcFilter> statusValues;

        if(fieldsIndexMap.containsKey("QcStatus"))
        {
            String qcStatusStr = values[fieldsIndexMap.get("QcStatus")];
            statusValues = qcFiltersFromString(qcStatusStr);
        }
        else
        {
            statusValues = calcQcStatus(totalFragments, duplicateFragments, splicedGeneCount);
        }

        return ImmutableRnaStatistics.builder()
                .qcStatus(statusValues)
                .totalFragments(totalFragments)
                .duplicateFragments(duplicateFragments)
                .splicedFragmentPerc(getDoubleValue(fieldsIndexMap, "SplicedFragmentPerc", values))
                .unsplicedFragmentPerc(getDoubleValue(fieldsIndexMap, "UnsplicedFragmentPerc", values))
                .altFragmentPerc(getDoubleValue(fieldsIndexMap, "AltFragmentPerc", values))
                .chimericFragmentPerc(getDoubleValue(fieldsIndexMap, "ChimericFragmentPerc", values))
                .splicedGeneCount(FileReaderUtils.getIntValue(fieldsIndexMap, "SplicedGeneCount", 0, values))
                .readLength(getIntValue(fieldsIndexMap, "ReadLength", values))
                .fragmentLength5thPercent(getDoubleValue(fieldsIndexMap, "FragLength5th", values))
                .fragmentLength50thPercent(getDoubleValue(fieldsIndexMap, "FragLength50th", values))
                .fragmentLength95thPercent(getDoubleValue(fieldsIndexMap, "FragLength95th", values))
                .enrichedGenePercent(getDoubleValue(fieldsIndexMap, "EnrichedGenePercent", values))
                .medianGCRatio(getDoubleValue(fieldsIndexMap, "MedianGCRatio", values))
                .forwardStrandPercent(getDoubleValue(fieldsIndexMap, "ForwardStrandPercent", values))
                .build();
    }

    public static List<RnaQcFilter> calcQcStatus(long totalFragments, long duplicateFragments, int splicedGenes)
    {
        return calcQcStatus(totalFragments, duplicateFragments, splicedGenes, LOW_COVERAGE_THRESHOLD, SPLICE_GENE_THRESHOLD);
    }

    public static List<RnaQcFilter> calcQcStatus(
            long totalFragments, long duplicateFragments, int splicedGenes, int lowCoverageThreshold, int splicedGeneThreshold)
    {
        if(totalFragments - duplicateFragments < lowCoverageThreshold)
            return List.of(FAIL_LOW_COVERAGE);

        List<RnaQcFilter> statusValues = Lists.newArrayListWithExpectedSize(2);

        double duplicatePerc = totalFragments > 0 ? duplicateFragments / (double)totalFragments : 0;

        if(duplicatePerc > HIGH_DUPLICATES_THRESHOLD)
            statusValues.add(WARN_DUPLICATE_RATE);

        if(splicedGenes >= 0 && splicedGenes < splicedGeneThreshold)
            statusValues.add(WARN_SPLICED_GENE_COVERAGE);

        return statusValues.isEmpty() ? List.of(PASS) : statusValues;
    }
}
