package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_ALT_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_CHIMERIC_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_DUPLICATE_FRAGS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_ENRICHED_GENE_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FORWARD_STRAND_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_50TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_5TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_95TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_MEDIAN_GC_RATIO;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_QC_STATUS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_READ_LENGTH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_SPLICED_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_UNSPLICED_FRAG_PERC;

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
        comparer = new IsofoxSummaryComparer(new ComparConfig());
        builder = TestIsofoxSummaryDataBuilder.BUILDER;
        IsofoxSummaryData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.ofEntries(
                Map.entry(FLD_QC_STATUS, b -> b.qcStatus = alternateValueSource.RnaStatistics().qcStatus()),
                Map.entry(FLD_TOTAL_FRAGS, b -> b.totalFragments = alternateValueSource.RnaStatistics().totalFragments()),
                Map.entry(FLD_DUPLICATE_FRAGS, b -> b.duplicateFragments = alternateValueSource.RnaStatistics().duplicateFragments()),
                Map.entry(FLD_SPLICED_FRAG_PERC, b -> b.splicedFragmentPerc = alternateValueSource.RnaStatistics().splicedFragmentPerc()),
                Map.entry(FLD_UNSPLICED_FRAG_PERC, b -> b.unsplicedFragmentPerc = alternateValueSource.RnaStatistics().unsplicedFragmentPerc()),
                Map.entry(FLD_ALT_FRAG_PERC, b -> b.altFragmentPerc = alternateValueSource.RnaStatistics().altFragmentPerc()),
                Map.entry(FLD_CHIMERIC_FRAG_PERC, b -> b.chimericFragmentPerc = alternateValueSource.RnaStatistics().chimericFragmentPerc()),
                Map.entry(FLD_READ_LENGTH, b -> b.readLength = alternateValueSource.RnaStatistics().readLength()),
                Map.entry(FLD_FRAG_LENGTH_5TH, b -> b.fragmentLength5thPercent = alternateValueSource.RnaStatistics().fragmentLength5thPercent()),
                Map.entry(FLD_FRAG_LENGTH_50TH, b -> b.fragmentLength50thPercent = alternateValueSource.RnaStatistics().fragmentLength50thPercent()),
                Map.entry(FLD_FRAG_LENGTH_95TH, b -> b.fragmentLength95thPercent = alternateValueSource.RnaStatistics().fragmentLength95thPercent()),
                Map.entry(FLD_ENRICHED_GENE_PERC, b -> b.enrichedGenePercent = alternateValueSource.RnaStatistics().enrichedGenePercent()),
                Map.entry(FLD_MEDIAN_GC_RATIO, b -> b.medianGCRatio = alternateValueSource.RnaStatistics().medianGCRatio()),
                Map.entry(FLD_FORWARD_STRAND_PERC, b -> b.forwardStrandPercent = alternateValueSource.RnaStatistics().forwardStrandPercent())
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
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