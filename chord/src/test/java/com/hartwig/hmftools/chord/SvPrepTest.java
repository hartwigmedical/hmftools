package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestUtils.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestUtils.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.TMP_OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.sv.SvPrep;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SvPrepTest
{
    @Before
    public void setup()
    {
        new File(TMP_OUTPUT_DIR).mkdir();
    }

    @After
    public void teardown() throws IOException
    {
        FileUtils.deleteDirectory(new File(TMP_OUTPUT_DIR));
    }

    @Test
    public void canPrepStructuralVariants()
    {
        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(List.of(MINIMAL_SAMPLE))
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .build();

        SvPrep prep = new SvPrep(config);

        List<MutContextCount> actualContextCounts = prep.countMutationContexts(MINIMAL_SAMPLE);
        //actualContextCounts.forEach(System.out::println);

        List<MutContextCount> expectedContextCounts = List.of(
                new MutContextCount("DEL_0e00_1e03_bp", 1),
                new MutContextCount("DEL_1e03_1e04_bp", 0),
                new MutContextCount("DEL_1e04_1e05_bp", 1),
                new MutContextCount("DEL_1e05_1e06_bp", 0),
                new MutContextCount("DEL_1e06_1e07_bp", 0),
                new MutContextCount("DEL_1e07_Inf_bp",  0),
                new MutContextCount("DUP_0e00_1e03_bp", 2),
                new MutContextCount("DUP_1e03_1e04_bp", 0),
                new MutContextCount("DUP_1e04_1e05_bp", 0),
                new MutContextCount("DUP_1e05_1e06_bp", 0),
                new MutContextCount("DUP_1e06_1e07_bp", 0),
                new MutContextCount("DUP_1e07_Inf_bp",  0),
                new MutContextCount("INV_0e00_1e03_bp", 0),
                new MutContextCount("INV_1e03_1e04_bp", 0),
                new MutContextCount("INV_1e04_1e05_bp", 0),
                new MutContextCount("INV_1e05_1e06_bp", 0),
                new MutContextCount("INV_1e06_1e07_bp", 0),
                new MutContextCount("INV_1e07_Inf_bp",  1),
                new MutContextCount("TRA",              1)
        );

        assertEquals(expectedContextCounts, actualContextCounts);
    }

    @Test
    public void canPrepStructuralVariantsFromEmptyVcf()
    {
        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(List.of(EMPTY_SAMPLE))
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .build();

        SvPrep prep = new SvPrep(config);

        List<MutContextCount> contextCounts = prep.countMutationContexts(EMPTY_SAMPLE);

        int contextCountTotal = 0;
        for(MutContextCount contextCount : contextCounts)
        {
            contextCountTotal += contextCount.mCount;
        }

        assertEquals(0, contextCountTotal);
    }
}
