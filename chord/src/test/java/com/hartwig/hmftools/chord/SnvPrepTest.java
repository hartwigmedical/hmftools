package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestDataPaths.DUMMY_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.NON_STANDARD_NUC_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.NON_STANDARD_NUC_GENOME_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.TMP_OUTPUT_DIR;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.snv.SnvPrep;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
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
        Configurator.setRootLevel(Level.DEBUG);

        ChordConfig config = new ChordConfig.Builder()
                .snvIndelVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SOMATIC_VCF_SUFFIX)
                .refGenomeFile(DUMMY_GENOME_FASTA)
                .outputDir(TMP_OUTPUT_DIR)
                .build();

        SnvPrep prep = new SnvPrep(config);

        List<MutContextCount> actualContextCounts = prep.countMutationContexts(MINIMAL_SAMPLE);
        //actualContextCounts.forEach(System.out::println);

        // There are 96 contexts in total, but only test the first few
        List<MutContextCount> firstExpectedContextCounts = List.of(
            new MutContextCount("A[C>A]A", 10),
            new MutContextCount("A[C>A]C", 6),
            new MutContextCount("A[C>A]G", 8),
            new MutContextCount("A[C>A]T", 3),
            new MutContextCount("C[C>A]A", 4),
            new MutContextCount("C[C>A]C", 4),
            new MutContextCount("C[C>A]G", 6),
            new MutContextCount("C[C>A]T", 5),
            new MutContextCount("G[C>A]A", 4),
            new MutContextCount("G[C>A]C", 4),
            new MutContextCount("G[C>A]G", 9),
            new MutContextCount("G[C>A]T", 5),
            new MutContextCount("T[C>A]A", 6),
            new MutContextCount("T[C>A]C", 4),
            new MutContextCount("T[C>A]G", 0),
            new MutContextCount("T[C>A]T", 6),
            new MutContextCount("A[C>G]A", 8),
            new MutContextCount("A[C>G]C", 1),
            new MutContextCount("A[C>G]G", 8),
            new MutContextCount("A[C>G]T", 11)
        );

        List<MutContextCount> firstActualContextCounts = new ArrayList<>();
        for(int i = 0; i < firstExpectedContextCounts.size(); i++)
        {
            firstActualContextCounts.add(actualContextCounts.get(i));
        }

        assertEquals(firstExpectedContextCounts, firstActualContextCounts);
    }

    @Test
    public void canPrepSnvsWithNonStandardNucleotides()
    {
        Configurator.setRootLevel(Level.DEBUG);

        ChordConfig config = new ChordConfig.Builder()
                .snvIndelVcfFile(INPUT_VCF_DIR + NON_STANDARD_NUC_GENOME_SAMPLE + PURPLE_SOMATIC_VCF_SUFFIX)
                .refGenomeFile(NON_STANDARD_NUC_GENOME_FASTA)
                .outputDir(TMP_OUTPUT_DIR)
                .build();

        SnvPrep prep = new SnvPrep(config);

        List<MutContextCount> actualContextCounts = prep.countMutationContexts(NON_STANDARD_NUC_GENOME_SAMPLE);

        assertTrue(actualContextCounts.stream().mapToInt(x -> x.mCount).sum() > 0);
    }

    @Test(expected = IllegalStateException.class)
    public void providingWrongVcfTypeThrowsError() throws NoSuchFileException
    {
        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(MINIMAL_SAMPLE)
                .snvIndelVcfFile(INPUT_VCF_DIR + MINIMAL_SAMPLE + PURPLE_SV_VCF_SUFFIX)
                .build();

        new SnvPrep(config).loadVariants(MINIMAL_SAMPLE);
    }
}
