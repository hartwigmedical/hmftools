package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
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

    public abstract String qcStatus();

    private static final int LOW_COVERAGE_THRESHOLD = 2500000;
    private static final int SPLICE_GENE_THRESHOLD = 17000;
    private static final double HIGH_DUPLICATES_THRESHOLD = 0.9;

    public static final String QC_PASS = "PASS";
    public static final String QC_FAIL_LOW_COVERAGE = "FAIL_LOW_COVERAGE";
    public static final String QC_WARN_DUPLICATES = "WARN_DUPLICATE_RATE";
    public static final String QC_WARN_LOW_SPLICE_GENES = "WARN_SPLICED_GENE_COVERAGE";

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
                .add(qcStatus())
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

        String qcStatus;

        if(fieldsIndexMap.containsKey("QcStatus"))
        {
            qcStatus = values[fieldsIndexMap.get("QcStatus")];
        }
        else
        {
            qcStatus = calcQcStatus(totalFragments, duplicateFragments, splicedGeneCount);
        }

        return ImmutableRnaStatistics.builder()
                .qcStatus(qcStatus)
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

    public static String calcQcStatus(long totalFragments, long duplicateFragments, int splicedGenes)
    {
        if(totalFragments - duplicateFragments < LOW_COVERAGE_THRESHOLD)
            return QC_FAIL_LOW_COVERAGE;

        List<String> warnings = Lists.newArrayList();

        double duplicatePerc = totalFragments > 0 ? duplicateFragments / (double)totalFragments : 0;

        if(duplicatePerc > HIGH_DUPLICATES_THRESHOLD)
            warnings.add(QC_WARN_DUPLICATES);

        if(splicedGenes >= 0 && splicedGenes < SPLICE_GENE_THRESHOLD)
            warnings.add(QC_WARN_LOW_SPLICE_GENES);

        if(warnings.isEmpty())
            return QC_PASS;

        StringJoiner sj = new StringJoiner(";");
        warnings.forEach(x -> sj.add(x));
        return sj.toString();
    }

}
