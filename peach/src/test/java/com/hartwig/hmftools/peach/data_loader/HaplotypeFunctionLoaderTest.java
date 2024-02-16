package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.peach.effect.HaplotypeFunction;
import com.hartwig.hmftools.peach.effect.ImmutableHaplotypeFunction;

import org.junit.Test;

public class HaplotypeFunctionLoaderTest
{
    @Test
    public void testLoad()
    {
        String filePath = getTestResourcePath("function.tsv");
        List<HaplotypeFunction> functions = HaplotypeFunctionLoader.loadFunctions(filePath);

        assertEquals(11, functions.size());

        HaplotypeFunction expectedFunction0 =
                ImmutableHaplotypeFunction.builder().geneName("DPYD").haplotypeName("*1").function("Normal Function").build();
        assertEquals(expectedFunction0, functions.get(0));

        HaplotypeFunction expectedFunction10 =
                ImmutableHaplotypeFunction.builder().geneName("UGT1A1").haplotypeName("*6").function("Reduced Function").build();
        assertEquals(expectedFunction10, functions.get(10));
    }
}
