package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestUtils.HUMAN_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestUtils.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestUtils.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.TMP_OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.indel.IndelPrep;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
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

    @Ignore
    @Test
    public void canPrepIndels()
    {
        Configurator.setRootLevel(Level.DEBUG);

        ChordConfig config = new ChordConfig.Builder()
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .refGenomeFile(HUMAN_GENOME_FASTA)
                .build();

        IndelPrep prep = new IndelPrep(config);
        List<MutContextCount> actualContextCounts = prep.countMutationContexts(MINIMAL_SAMPLE);
        //actualContextCounts.forEach(System.out::println);

        List<MutContextCount> expectedContextCounts = List.of(
                new MutContextCount("del.rep.len.1", 7),
                new MutContextCount("del.rep.len.2", 1),
                new MutContextCount("del.rep.len.3", 0),
                new MutContextCount("del.rep.len.4", 1),
                new MutContextCount("del.rep.len.5", 4),
                new MutContextCount("ins.rep.len.1", 4),
                new MutContextCount("ins.rep.len.2", 2),
                new MutContextCount("ins.rep.len.3", 1),
                new MutContextCount("ins.rep.len.4", 0),
                new MutContextCount("ins.rep.len.5", 1),
                new MutContextCount("del.mh.bimh.1", 0),
                new MutContextCount("del.mh.bimh.2", 1),
                new MutContextCount("del.mh.bimh.3", 1),
                new MutContextCount("del.mh.bimh.4", 0),
                new MutContextCount("del.mh.bimh.5", 1),
                new MutContextCount("ins.mh.bimh.1", 0),
                new MutContextCount("ins.mh.bimh.2", 0),
                new MutContextCount("ins.mh.bimh.3", 2),
                new MutContextCount("ins.mh.bimh.4", 1),
                new MutContextCount("ins.mh.bimh.5", 1),
                new MutContextCount("del.none.len.1", 4),
                new MutContextCount("del.none.len.2", 0),
                new MutContextCount("del.none.len.3", 0),
                new MutContextCount("del.none.len.4", 0),
                new MutContextCount("del.none.len.5", 0),
                new MutContextCount("ins.none.len.1", 1),
                new MutContextCount("ins.none.len.2", 0),
                new MutContextCount("ins.none.len.3", 0),
                new MutContextCount("ins.none.len.4", 0),
                new MutContextCount("ins.none.len.5", 1)
        );

        assertEquals(expectedContextCounts, actualContextCounts);
    }
}
