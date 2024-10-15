package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestUtils.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.HUMAN_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestUtils.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestUtils.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.TMP_OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.snv.SnvPrep;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SnvPrepTest
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
    public void canPrepSnvs()
    {
        ChordConfig config = new ChordConfig.Builder()
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .refGenomeFile(HUMAN_GENOME_FASTA)
                .build();

        SnvPrep prep = new SnvPrep(config);

        List<MutContextCount> actualContextCounts = prep.countMutationContexts(MINIMAL_SAMPLE);
        //actualContextCounts.forEach(System.out::println);

        // There are 96 contexts in total, but only test the first few
        List<MutContextCount> firstExpectedContextCounts = List.of(
                new MutContextCount("A[C>A]A", 8),
                new MutContextCount("A[C>A]C", 1),
                new MutContextCount("A[C>A]G", 3),
                new MutContextCount("A[C>A]T", 2),
                new MutContextCount("C[C>A]A", 87),
                new MutContextCount("C[C>A]C", 10),
                new MutContextCount("C[C>A]G", 21),
                new MutContextCount("C[C>A]T", 48),
                new MutContextCount("G[C>A]A", 3),
                new MutContextCount("G[C>A]C", 1),
                new MutContextCount("G[C>A]G", 0),
                new MutContextCount("G[C>A]T", 1),
                new MutContextCount("T[C>A]A", 15),
                new MutContextCount("T[C>A]C", 7),
                new MutContextCount("T[C>A]G", 1),
                new MutContextCount("T[C>A]T", 12),
                new MutContextCount("A[C>G]A", 2),
                new MutContextCount("A[C>G]C", 4),
                new MutContextCount("A[C>G]G", 0),
                new MutContextCount("A[C>G]T", 3)
        );

        List<MutContextCount> firstActualContextCounts = new ArrayList<>();
        for(int i = 0; i < firstExpectedContextCounts.size(); i++)
        {
            firstActualContextCounts.add(actualContextCounts.get(i));
        }

        assertEquals(firstExpectedContextCounts, firstActualContextCounts);
    }
}
