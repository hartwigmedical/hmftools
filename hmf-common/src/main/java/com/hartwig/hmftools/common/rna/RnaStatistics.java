package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.DELIMITER;

import java.util.StringJoiner;

import org.immutables.value.Value;

@Value.Immutable
public abstract class RnaStatistics
{
    public abstract int totalFragments();
    public abstract int duplicateFragments();
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

    public static final int LOW_COVERAGE_THRESHOLD = 2500000;

    public String qcStatus()
    {
        return totalFragments() - duplicateFragments() >= LOW_COVERAGE_THRESHOLD ? "PASS" : "FAIL_LOW_COVERAGE";
    }

    public static String csvHeader()
    {
        return "SampleId,TotalFragments,DuplicateFragments"
                + ",SplicedFragmentPerc,UnsplicedFragmentPerc,AltFragmentPerc,ChimericFragmentPerc"
                + ",ReadLength,FragLength5th,FragLength50th,FragLength95th"
                + ",EnrichedGenePercent,MedianGCRatio";
    }

    public String toCsv(final String sampleId)
    {
        return new StringJoiner(DELIMITER)
                .add(sampleId)
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

    public static RnaStatistics fromCsv(final String input)
    {
        final String[] items = input.split(DELIMITER);
        int index = 1;

        return ImmutableRnaStatistics.builder()
                .totalFragments(Integer.parseInt(items[index++]))
                .duplicateFragments(Integer.parseInt(items[index++]))
                .splicedFragmentPerc(Double.parseDouble(items[index++]))
                .unsplicedFragmentPerc(Double.parseDouble(items[index++]))
                .altFragmentPerc(Double.parseDouble(items[index++]))
                .chimericFragmentPerc(Double.parseDouble(items[index++]))
                .readLength(Integer.parseInt(items[index++]))
                .fragmentLength5thPercent(Double.parseDouble(items[index++]))
                .fragmentLength50thPercent(Double.parseDouble(items[index++]))
                .fragmentLength95thPercent(Double.parseDouble(items[index++]))
                .enrichedGenePercent(Double.parseDouble(items[index++]))
                .medianGCRatio(Double.parseDouble(items[index++]))
                .build();
    }
}
