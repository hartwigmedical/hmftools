package com.hartwig.hmftools.viridian.app;

import static com.hartwig.hmftools.viridian.app.ViridianConfig.parseExtraVirusDecoyContigs;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.Set;

import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.profiling.extract.CandidateReadSource;

import org.junit.Test;

public class ViridianConfigTest
{
    @Test
    public void testParseExtraVirusDecoyContigsUnset()
    {
        assertEquals(Set.of(), parseExtraVirusDecoyContigs(null));
    }

    @Test
    public void testParseExtraVirusDecoyContigs()
    {
        assertEquals(Set.of("chrHPV16"), parseExtraVirusDecoyContigs("chrHPV16"));
        assertEquals(Set.of("chrHPV16", "chrHBV"), parseExtraVirusDecoyContigs("chrHPV16;chrHBV"));
    }

    @Test
    public void testParseExtraVirusDecoyContigsBlankContig()
    {
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs(""));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16;;chrHBV"));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16;"));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16; "));
    }

    @Test
    public void testConstructorCounts()
    {
        config(1, 1);
        assertThrows(UserInputError.class, () -> config(0, 1));
        assertThrows(UserInputError.class, () -> config(1, 0));
    }

    private static ViridianConfig config(int alignmentBatchSize, int threads)
    {
        return new ViridianConfig(
                "sample", "tumor.bam", "ref.fasta", "virus.fasta", "virus_genomes.tsv", "virus.fasta.img",
                "virus_groups.tsv", null, "esvee.vcf.gz", "purple", "bam_metrics", alignmentBatchSize, threads, "output/", null,
                CandidateReadSource.FULLY_UNMAPPED_AND_DECOY, Set.of(), false);
    }
}
