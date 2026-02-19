package com.hartwig.hmftools.compar.isofox;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.rna.ImmutableRnaStatistics;
import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestIsofoxSummaryDataBuilder
{
    public List<RnaQcFilter> qcStatus = List.of(RnaQcFilter.PASS);
    public long totalFragments = 1000000;
    public long duplicateFragments = 500000;
    public double splicedFragmentPerc = 0.7;
    public double unsplicedFragmentPerc = 0.3;
    public double altFragmentPerc = 0.1;
    public double chimericFragmentPerc = 0.1;
    public int readLength = 151;
    public double fragmentLength5thPercent = 31;
    public double fragmentLength50thPercent = 100;
    public double fragmentLength95thPercent = 150;
    public double enrichedGenePercent = 0.2;
    public double medianGCRatio = 0.5;
    public double forwardStrandPercent = 0.9;

    private static final Consumer<TestIsofoxSummaryDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.qcStatus = List.of(RnaQcFilter.FAIL_LOW_COVERAGE);
        b.totalFragments = 2000000;
        b.duplicateFragments = 1000000;
        b.splicedFragmentPerc = 0.6;
        b.unsplicedFragmentPerc = 0.4;
        b.altFragmentPerc = 0.2;
        b.chimericFragmentPerc = 0.2;
        b.readLength = 76;
        b.fragmentLength5thPercent = 25;
        b.fragmentLength50thPercent = 80;
        b.fragmentLength95thPercent = 120;
        b.enrichedGenePercent = 0.4;
        b.medianGCRatio = 0.6;
        b.forwardStrandPercent = 0.5;
    };

    public static final TestComparableItemBuilder<TestIsofoxSummaryDataBuilder, IsofoxSummaryData> BUILDER =
            new TestComparableItemBuilder<>(TestIsofoxSummaryDataBuilder::new, TestIsofoxSummaryDataBuilder::build, ALTERNATE_INITIALIZER);

    private IsofoxSummaryData build()
    {
        final RnaStatistics rnaStatistics = ImmutableRnaStatistics.builder()
                .qcStatus(qcStatus)
                .totalFragments(totalFragments)
                .duplicateFragments(duplicateFragments)
                .splicedFragmentPerc(splicedFragmentPerc)
                .unsplicedFragmentPerc(unsplicedFragmentPerc)
                .altFragmentPerc(altFragmentPerc)
                .chimericFragmentPerc(chimericFragmentPerc)
                .splicedGeneCount(-1)
                .readLength(readLength)
                .fragmentLength5thPercent(fragmentLength5thPercent)
                .fragmentLength50thPercent(fragmentLength50thPercent)
                .fragmentLength95thPercent(fragmentLength95thPercent)
                .enrichedGenePercent(enrichedGenePercent)
                .medianGCRatio(medianGCRatio)
                .forwardStrandPercent(forwardStrandPercent)
                .build();
        return new IsofoxSummaryData(rnaStatistics);
    }
}
