package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class RnaGeneDataTest extends ComparableItemTest<RnaGeneData, RnaGeneDataComparer, TestIsofoxGeneDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaGeneDataComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestIsofoxGeneDataBuilder.BUILDER;
        RnaGeneData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaGeneDataComparer.Fields.SplicedFragments.toString(), b -> b.splicedFragments = alternateValueSource.Expression.splicedFragments(),
                RnaGeneDataComparer.Fields.UnsplicedFragments.toString(), b -> b.unsplicedFragments = alternateValueSource.Expression.unsplicedFragments(),
                RnaGeneDataComparer.Fields.AdjTPM.toString(), b -> b.tpm = alternateValueSource.Expression.tpm()
        );
        nameToAlternateIndexInitializer = Map.of(
                FLD_GENE_NAME, b -> b.geneName = alternateValueSource.Expression.geneName()
        );
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
