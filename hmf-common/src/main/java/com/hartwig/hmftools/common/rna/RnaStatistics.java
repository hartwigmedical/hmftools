package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

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

    public abstract int readLength();

    // 5th, 50th and 95th percentile intronic fragment length
    public abstract double fragmentLength5thPercent();
    public abstract double fragmentLength50thPercent();
    public abstract double fragmentLength95thPercent();

    // proportion of fragments in 7 highly expressed genes
    public abstract double enrichedGenePercent();

    // Median GC (excluding 7 highly expressed genes)
    public abstract double medianGCRatio();

    public abstract String qcStatus();

    private static final int LOW_COVERAGE_THRESHOLD = 2500000;
    private static final int SPLICE_GENE_THRESHOLD = 17000;
    private static final double HIGH_DUPLICATES_THRESHOLD = 0.9;

    public static final String QC_PASS = "PASS";
    public static final String QC_FAIL_LOW_COVERAGE = "FAIL_LOW_COVERAGE";
    public static final String QC_WARN_DUPLICATES = "WARN_DUPLICATE_RATE";
    public static final String QC_WARN_LOW_SPLICE_GENES = "WARN_SPLICED_GENE_COVERAGE";

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
                .add("ReadLength")
                .add("FragLength5th")
                .add("FragLength50th")
                .add("FragLength95th")
                .add("EnrichedGenePercent")
                .add("MedianGCRatio").toString();
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
                .add(String.valueOf(readLength()))
                .add(String.format("%.0f", fragmentLength5thPercent()))
                .add(String.format("%.0f", fragmentLength50thPercent()))
                .add(String.format("%.0f", fragmentLength95thPercent()))
                .add(String.format("%.3f", enrichedGenePercent()))
                .add(String.format("%.3f", medianGCRatio()))
                .toString();
    }

    public static RnaStatistics fromCsv(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);

        final String[] values = lines.get(1).split(DELIMITER);

        ImmutableRnaStatistics.Builder builder = ImmutableRnaStatistics.builder();

        long totalFragments = Long.parseLong(values[fieldsIndexMap.get("TotalFragments")]);
        long duplicateFragments = Long.parseLong(values[fieldsIndexMap.get("DuplicateFragments")]);

        String qcStatus;

        if(fieldsIndexMap.containsKey("QcStatus"))
        {
            qcStatus = values[fieldsIndexMap.get("QcStatus")];
        }
        else
        {
            qcStatus = calcQcStatus(totalFragments, duplicateFragments, -1);
        }

        return ImmutableRnaStatistics.builder()
                .qcStatus(qcStatus)
                .totalFragments(totalFragments)
                .duplicateFragments(duplicateFragments)
                .splicedFragmentPerc(Double.parseDouble(values[fieldsIndexMap.get("SplicedFragmentPerc")]))
                .unsplicedFragmentPerc(Double.parseDouble(values[fieldsIndexMap.get("UnsplicedFragmentPerc")]))
                .altFragmentPerc(Double.parseDouble(values[fieldsIndexMap.get("AltFragmentPerc")]))
                .chimericFragmentPerc(Double.parseDouble(values[fieldsIndexMap.get("ChimericFragmentPerc")]))
                .readLength(Integer.parseInt(values[fieldsIndexMap.get("ReadLength")]))
                .fragmentLength5thPercent(Double.parseDouble(values[fieldsIndexMap.get("FragLength5th")]))
                .fragmentLength50thPercent(Double.parseDouble(values[fieldsIndexMap.get("FragLength50th")]))
                .fragmentLength95thPercent(Double.parseDouble(values[fieldsIndexMap.get("FragLength95th")]))
                .enrichedGenePercent(Double.parseDouble(values[fieldsIndexMap.get("EnrichedGenePercent")]))
                .medianGCRatio(Double.parseDouble(values[fieldsIndexMap.get("MedianGCRatio")]))
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
