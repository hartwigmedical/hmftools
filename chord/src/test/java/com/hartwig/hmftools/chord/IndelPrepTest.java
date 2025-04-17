package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestDataPaths.DUMMY_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.TMP_OUTPUT_DIR;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.indel.IndelPrep;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class IndelPrepTest
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
    public void canPrepIndels()
    {
        Configurator.setRootLevel(Level.DEBUG);

        ChordConfig config = new ChordConfig.Builder()
                .snvIndelVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SOMATIC_VCF_SUFFIX)
                .refGenomeFile(DUMMY_GENOME_FASTA)
                .outputDir(TMP_OUTPUT_DIR)
                .build();

        IndelPrep prep = new IndelPrep(config);
        List<MutContextCount> actualContextCounts = prep.countMutationContexts(MINIMAL_SAMPLE);
        //actualContextCounts.forEach(System.out::println);

        List<MutContextCount> expectedContextCounts = List.of(
            new MutContextCount("del.rep.len.1", 30),
            new MutContextCount("del.rep.len.2", 3),
            new MutContextCount("del.rep.len.3", 1),
            new MutContextCount("del.rep.len.4", 1),
            new MutContextCount("del.rep.len.5", 0),
            new MutContextCount("ins.rep.len.1", 25),
            new MutContextCount("ins.rep.len.2", 2),
            new MutContextCount("ins.rep.len.3", 1),
            new MutContextCount("ins.rep.len.4", 1),
            new MutContextCount("ins.rep.len.5", 0),
            new MutContextCount("del.mh.bimh.1", 8),
            new MutContextCount("del.mh.bimh.2", 11),
            new MutContextCount("del.mh.bimh.3", 0),
            new MutContextCount("del.mh.bimh.4", 0),
            new MutContextCount("del.mh.bimh.5", 1),
            new MutContextCount("ins.mh.bimh.1", 9),
            new MutContextCount("ins.mh.bimh.2", 12),
            new MutContextCount("ins.mh.bimh.3", 2),
            new MutContextCount("ins.mh.bimh.4", 0),
            new MutContextCount("ins.mh.bimh.5", 0),
            new MutContextCount("del.none.len.1", 75),
            new MutContextCount("del.none.len.2", 43),
            new MutContextCount("del.none.len.3", 26),
            new MutContextCount("del.none.len.4", 7),
            new MutContextCount("del.none.len.5", 8),
            new MutContextCount("ins.none.len.1", 99),
            new MutContextCount("ins.none.len.2", 55),
            new MutContextCount("ins.none.len.3", 35),
            new MutContextCount("ins.none.len.4", 14),
            new MutContextCount("ins.none.len.5", 13)
        );

        assertEquals(expectedContextCounts, actualContextCounts);
    }

    @Test(expected = IllegalStateException.class)
    public void providingWrongVcfTypeThrowsError() throws NoSuchFileException
    {
        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(MINIMAL_SAMPLE)
                .snvIndelVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SOMATIC_VCF_SUFFIX)
                .build();

        new IndelPrep(config).loadVariants(MINIMAL_SAMPLE);
    }
}
