package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class RnaNovelSpliceJunctionDataTest
        extends ComparableItemTest<RnaNovelSpliceJunctionData, RnaNovelSpliceJunctionComparer, TestNovelSpliceJunctionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaNovelSpliceJunctionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestNovelSpliceJunctionDataBuilder.BUILDER;
        RnaNovelSpliceJunctionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaNovelSpliceJunctionComparer.Fields.Type.toString(), b -> b.type = alternateValueSource.Type,
                RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString(), b -> b.fragmentCount = alternateValueSource.FragmentCount,
                RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString(), b -> b.regionStart = alternateValueSource.RegionTypeStart,
                RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString(), b -> b.regionEnd = alternateValueSource.RegionTypeEnd);

        nameToAlternateIndexInitializer = Map.of(
                FLD_GENE_NAME, b -> b.geneName = alternateValueSource.GeneName,
                FLD_CHROMOSOME, b -> {
                    b.chromosome = alternateValueSource.Chromosome;
                },
                FLD_ALT_SJ_POS_START, b -> {
                    b.junctionStart = alternateValueSource.JunctionStart;
                },
                FLD_ALT_SJ_POS_END, b -> {
                    b.junctionEnd = alternateValueSource.JunctionEnd;
                });
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
