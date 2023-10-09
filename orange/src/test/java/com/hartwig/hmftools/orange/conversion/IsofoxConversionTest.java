package com.hartwig.hmftools.orange.conversion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;

import org.junit.Test;

public class IsofoxConversionTest
{
    @Test
    public void shouldConvertCancerFlagValuesToNull()
    {
        ImmutableGeneExpression expressionGood =
                IsofoxTestFactory.geneExpressionBuilder().medianTpmCancer(10D).percentileCancer(20D).build();
        final GeneExpression convertedGood = IsofoxConversion.convert(expressionGood);
        assertNotNull(convertedGood.medianTpmCancer());
        assertNotNull(convertedGood.percentileCancer());

        assertEquals(10D, convertedGood.medianTpmCancer(), 0D);
        assertEquals(20D, convertedGood.percentileCancer(), 0D);

        ImmutableGeneExpression expressionFlagged =
                IsofoxTestFactory.geneExpressionBuilder().medianTpmCancer(-1D).percentileCancer(-1D).build();
        final GeneExpression convertedFlagged = IsofoxConversion.convert(expressionFlagged);
        assertNull(convertedFlagged.medianTpmCancer());
        assertNull(convertedFlagged.percentileCancer());
    }
}
