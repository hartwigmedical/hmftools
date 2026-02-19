package com.hartwig.hmftools.orange.conversion;

import static com.hartwig.hmftools.orange.conversion.IsofoxConversion.convert;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.orange.algo.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;

import org.junit.Test;

public class IsofoxConversionTest
{
    @Test
    public void shouldConvertCancerFlagValuesToNull()
    {
        ImmutableGeneExpression expressionGood =
                IsofoxTestFactory.geneExpressionBuilder().medianTpmCancer(10D).percentileCancer(20D).build();
        final GeneExpression convertedGood = convert(expressionGood);
        assertNotNull(convertedGood.medianTpmCancer());
        assertNotNull(convertedGood.percentileCancer());

        assertEquals(10D, convertedGood.medianTpmCancer(), 0D);
        assertEquals(20D, convertedGood.percentileCancer(), 0D);

        ImmutableGeneExpression expressionFlagged =
                IsofoxTestFactory.geneExpressionBuilder().medianTpmCancer(-1D).percentileCancer(-1D).build();
        final GeneExpression convertedFlagged = convert(expressionFlagged);
        assertNull(convertedFlagged.medianTpmCancer());
        assertNull(convertedFlagged.percentileCancer());
    }

    @Test
    public void shouldConvertFusionsWithoutGeneNames()
    {
        ImmutableRnaFusion fusion = IsofoxTestFactory.rnaFusionBuilder()
                .chromosomeUp("1")
                .chromosomeDown("1")
                .positionUp(100)
                .positionDown(200)
                .junctionTypeUp("UNKNOWN")
                .junctionTypeDown("UNKNOWN")
                .svType(StructuralVariantType.DEL)
                .build();

        RnaFusion convertedFusion = convert(fusion);
        assertNull(convertedFusion.geneStart());
        assertNull(convertedFusion.geneEnd());
        assertEquals("::", convertedFusion.display());
    }
}
