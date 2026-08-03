package com.hartwig.hmftools.compar.isofox;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class IsofoxSummaryDataTest extends ComparableItemTest<IsofoxSummaryData, IsofoxSummaryComparer, TestIsofoxSummaryDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new IsofoxSummaryComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestIsofoxSummaryDataBuilder.BUILDER;
        IsofoxSummaryData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.ofEntries(
                Map.entry(IsofoxSummaryComparer.Fields.QcStatus.toString(), b -> b.qcStatus = alternateValueSource.Statistics.qcStatus()),
                Map.entry(IsofoxSummaryComparer.Fields.TotalFragments.toString(), b -> b.totalFragments = alternateValueSource.Statistics.totalFragments()),
                Map.entry(IsofoxSummaryComparer.Fields.DuplicateFragments.toString(), b -> b.duplicateFragments = alternateValueSource.Statistics.duplicateFragments()),
                Map.entry(IsofoxSummaryComparer.Fields.SplicedFragmentPerc.toString(), b -> b.splicedFragmentPerc = alternateValueSource.Statistics.splicedFragmentPerc()),
                Map.entry(IsofoxSummaryComparer.Fields.UnsplicedFragmentPerc.toString(), b -> b.unsplicedFragmentPerc = alternateValueSource.Statistics.unsplicedFragmentPerc()),
                Map.entry(IsofoxSummaryComparer.Fields.AltFragmentPerc.toString(), b -> b.altFragmentPerc = alternateValueSource.Statistics.altFragmentPerc()),
                Map.entry(IsofoxSummaryComparer.Fields.ChimericFragmentPerc.toString(), b -> b.chimericFragmentPerc = alternateValueSource.Statistics.chimericFragmentPerc()),
                Map.entry(IsofoxSummaryComparer.Fields.FragLength5th.toString(), b -> b.fragmentLength5thPercent = alternateValueSource.Statistics.fragmentLength5thPercent()),
                Map.entry(IsofoxSummaryComparer.Fields.FragLength50th.toString(), b -> b.fragmentLength50thPercent = alternateValueSource.Statistics.fragmentLength50thPercent()),
                Map.entry(IsofoxSummaryComparer.Fields.FragLength95th.toString(), b -> b.fragmentLength95thPercent = alternateValueSource.Statistics.fragmentLength95thPercent()),
                Map.entry(IsofoxSummaryComparer.Fields.MedianGCRatio.toString(), b -> b.medianGCRatio = alternateValueSource.Statistics.medianGCRatio()),
                Map.entry(IsofoxSummaryComparer.Fields.ForwardStrandPercent.toString(), b -> b.forwardStrandPercent = alternateValueSource.Statistics.forwardStrandPercent())
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }
}