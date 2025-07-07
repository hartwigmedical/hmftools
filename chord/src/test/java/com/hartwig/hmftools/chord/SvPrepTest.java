package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestDataPaths.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.TMP_OUTPUT_DIR;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.sv.SvPrep;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
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
        Configurator.setRootLevel(Level.DEBUG);

        ChordConfig config = new ChordConfig.Builder()
                .svVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SV_VCF_SUFFIX)
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
                .sampleIds(EMPTY_SAMPLE)
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

    @Test(expected = IllegalStateException.class)
    public void providingWrongVcfTypeThrowsError() throws NoSuchFileException
    {
        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(MINIMAL_SAMPLE)
                .svVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SOMATIC_VCF_SUFFIX)
                .build();

        new SvPrep(config).loadVariants(MINIMAL_SAMPLE);
    }
}
