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

        dup("PIK3R1", "D464_Y470dup");
        variant = duplication.calculateVariant(genome);
        checkSingleHotspot(variant, "T", "TGATAGATTATATGAAGAATAT", "chr5", 68293798);
    }

    private Duplication dup(String gene, String variant)
    {
        return transval.variationParser().parseDuplication(gene,variant);
    }
}
