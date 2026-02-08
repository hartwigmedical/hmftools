package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.FAIL_LOW_COVERAGE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.PASS;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.WARN_DUPLICATE_RATE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.WARN_SPLICED_GENE_COVERAGE;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersFromString;
import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersToString;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
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

public final class RnaStatisticFile
{
    public static final int LOW_COVERAGE_THRESHOLD = 2500000;
    public static final int SPLICE_GENE_THRESHOLD = 17000;
    private static final double HIGH_DUPLICATES_THRESHOLD = 0.9;

    public static final String SUMMARY_FILE_ID = "summary.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + SUMMARY_FILE_ID;
    }

    public enum Column
    {
        SampleId,
        QcStatus,
        TotalFragments,
        DuplicateFragments,
        SplicedFragmentPerc,
        UnsplicedFragmentPerc,
        AltFragmentPerc,
        ChimericFragmentPerc,
        SplicedGeneCount,
        ReadLength,
        FragLength5th,
        FragLength50th,
        FragLength95th,
        EnrichedGenePercent,
        MedianGCRatio,
        ForwardStrandPercent;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        for(Column column : RnaStatisticFile.Column.values())
        {
            sj.add(column.toString());
        }

        return sj.toString();
    }

    public static String writeLine(final String sampleId, final RnaStatistics statistics)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(sampleId);
        sj.add(qcFiltersToString(statistics.qcStatus()));
        sj.add(String.valueOf(statistics.totalFragments()));
        sj.add(String.valueOf(statistics.duplicateFragments()));
        sj.add(String.format("%.3f", statistics.splicedFragmentPerc()));
        sj.add(String.format("%.3f", statistics.unsplicedFragmentPerc()));
        sj.add(String.format("%.3f", statistics.altFragmentPerc()));
        sj.add(String.format("%.3f", statistics.chimericFragmentPerc()));
        sj.add(String.valueOf(statistics.splicedGeneCount()));
        sj.add(String.valueOf(statistics.readLength()));
        sj.add(String.format("%.0f", statistics.fragmentLength5thPercent()));
        sj.add(String.format("%.0f", statistics.fragmentLength50thPercent()));
        sj.add(String.format("%.0f", statistics.fragmentLength95thPercent()));
        sj.add(String.format("%.3f", statistics.enrichedGenePercent()));
        sj.add(String.format("%.3f", statistics.medianGCRatio()));
        sj.add(String.format("%.3f", statistics.forwardStrandPercent()));
        return sj.toString();
    }

    public static RnaStatistics fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        String fileDelim = inferHeaderDelimiter(header);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

        String[] values = lines.get(1).split(fileDelim);

        long totalFragments = getLongValue(fieldsIndexMap, Column.TotalFragments.toString(), values);
        long duplicateFragments = getLongValue(fieldsIndexMap, Column.DuplicateFragments.toString(), values);

        int splicedGeneCount = getIntValue(fieldsIndexMap, Column.SplicedGeneCount.toString(), values);

        List<RnaQcFilter> statusValues;

        if(fieldsIndexMap.containsKey(Column.QcStatus.toString()))
        {
            String qcStatusStr = values[fieldsIndexMap.get(Column.QcStatus.toString())];
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
                .splicedFragmentPerc(getDoubleValue(fieldsIndexMap, Column.SplicedFragmentPerc.toString(), values))
                .unsplicedFragmentPerc(getDoubleValue(fieldsIndexMap, Column.UnsplicedFragmentPerc.toString(), values))
                .altFragmentPerc(getDoubleValue(fieldsIndexMap, Column.AltFragmentPerc.toString(), values))
                .chimericFragmentPerc(getDoubleValue(fieldsIndexMap, Column.ChimericFragmentPerc.toString(), values))
                .splicedGeneCount(FileReaderUtils.getIntValue(fieldsIndexMap, Column.SplicedGeneCount.toString(), 0, values))
                .readLength(getIntValue(fieldsIndexMap, Column.ReadLength.toString(), values))
                .fragmentLength5thPercent(getDoubleValue(fieldsIndexMap, Column.FragLength5th.toString(), values))
                .fragmentLength50thPercent(getDoubleValue(fieldsIndexMap, Column.FragLength50th.toString(), values))
                .fragmentLength95thPercent(getDoubleValue(fieldsIndexMap, Column.FragLength95th.toString(), values))
                .enrichedGenePercent(getDoubleValue(fieldsIndexMap, Column.EnrichedGenePercent.toString(), values))
                .medianGCRatio(getDoubleValue(fieldsIndexMap, Column.MedianGCRatio.toString(), values))
                .forwardStrandPercent(getDoubleValue(fieldsIndexMap, Column.ForwardStrandPercent.toString(), values))
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
