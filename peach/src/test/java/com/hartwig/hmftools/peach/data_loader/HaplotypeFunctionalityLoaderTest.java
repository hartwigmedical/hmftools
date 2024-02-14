package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.peach.effect.HaplotypeFunctionality;
import com.hartwig.hmftools.peach.effect.ImmutableHaplotypeFunctionality;

import org.junit.Test;

public class HaplotypeFunctionalityLoaderTest
{
    @Test
    public void testLoad()
    {
        String filePath = getTestResourcePath("functionality.tsv");
        List<HaplotypeFunctionality> functionalities = HaplotypeFunctionalityLoader.loadFunctionalities(filePath);

        assertEquals(11, functionalities.size());

        HaplotypeFunctionality expectedFunctionality0 =
                ImmutableHaplotypeFunctionality.builder().geneName("DPYD").haplotypeName("*1").functionality("Normal Function").build();
        assertEquals(expectedFunctionality0, functionalities.get(0));

        HaplotypeFunctionality expectedFunctionality10 =
                ImmutableHaplotypeFunctionality.builder().geneName("UGT1A1").haplotypeName("*6").functionality("Reduced Function").build();
        assertEquals(expectedFunctionality10, functionalities.get(10));
    }
}
