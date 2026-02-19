package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class IsofoxGeneDataTest extends ComparableItemTest<IsofoxGeneData, IsofoxGeneDataComparer, TestIsofoxGeneDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new IsofoxGeneDataComparer(new ComparConfig());
        builder = TestIsofoxGeneDataBuilder.BUILDER;
        IsofoxGeneData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_SPLICED_FRAGS, b -> b.splicedFragments = alternateValueSource.GeneExpression().splicedFragments(),
                FLD_UNSPLICED_FRAGS, b -> b.unsplicedFragments = alternateValueSource.GeneExpression().unsplicedFragments(),
                FLD_ADJ_TPM, b -> b.tpm = alternateValueSource.GeneExpression().tpm()
        );
        nameToAlternateIndexInitializer = Map.of(
                FLD_GENE_NAME, b -> b.geneName = alternateValueSource.GeneExpression().geneName()
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
