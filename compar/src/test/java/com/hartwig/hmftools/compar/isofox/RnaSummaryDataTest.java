package com.hartwig.hmftools.compar.isofox;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class RnaSummaryDataTest extends ComparableItemTest<RnaSummaryData, RnaSummaryComparer, TestRnaSummaryDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaSummaryComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestRnaSummaryDataBuilder.BUILDER;
        RnaSummaryData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.ofEntries(
                Map.entry(RnaSummaryComparer.Fields.QcStatus.toString(), b -> b.qcStatus = alternateValueSource.Statistics.qcStatus()),
                Map.entry(RnaSummaryComparer.Fields.TotalFragments.toString(), b -> b.totalFragments = alternateValueSource.Statistics.totalFragments()),
                Map.entry(RnaSummaryComparer.Fields.DuplicateFragments.toString(), b -> b.duplicateFragments = alternateValueSource.Statistics.duplicateFragments()),
                Map.entry(RnaSummaryComparer.Fields.SplicedFragmentPerc.toString(), b -> b.splicedFragmentPerc = alternateValueSource.Statistics.splicedFragmentPerc()),
                Map.entry(RnaSummaryComparer.Fields.UnsplicedFragmentPerc.toString(), b -> b.unsplicedFragmentPerc = alternateValueSource.Statistics.unsplicedFragmentPerc()),
                Map.entry(RnaSummaryComparer.Fields.AltFragmentPerc.toString(), b -> b.altFragmentPerc = alternateValueSource.Statistics.altFragmentPerc()),
                Map.entry(RnaSummaryComparer.Fields.ChimericFragmentPerc.toString(), b -> b.chimericFragmentPerc = alternateValueSource.Statistics.chimericFragmentPerc()),
                Map.entry(RnaSummaryComparer.Fields.FragLength5th.toString(), b -> b.fragmentLength5thPercent = alternateValueSource.Statistics.fragmentLength5thPercent()),
                Map.entry(RnaSummaryComparer.Fields.FragLength50th.toString(), b -> b.fragmentLength50thPercent = alternateValueSource.Statistics.fragmentLength50thPercent()),
                Map.entry(RnaSummaryComparer.Fields.FragLength95th.toString(), b -> b.fragmentLength95thPercent = alternateValueSource.Statistics.fragmentLength95thPercent()),
                Map.entry(RnaSummaryComparer.Fields.MedianGCRatio.toString(), b -> b.medianGCRatio = alternateValueSource.Statistics.medianGCRatio()),
                Map.entry(RnaSummaryComparer.Fields.ForwardStrandPercent.toString(), b -> b.forwardStrandPercent = alternateValueSource.Statistics.forwardStrandPercent())
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