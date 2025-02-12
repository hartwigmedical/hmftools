package com.hartwig.hmftools.pave.transval;

import org.junit.Test;

public class DuplicationTest extends TransvalTestBase
{
    @Test
    public void calculateVariantTest()
    {
        Duplication duplication = dup("PIK3R1", "Y452dup");
        TransvalVariant variant = duplication.calculateVariant(genome);
        checkSingleHotspot(variant, "A", "ATAT", "chr5", 68293762);
    }

    private Duplication dup(String gene, String variant)
    {
        return transval.variationParser().parseDuplication(gene, variant);
    }
}
