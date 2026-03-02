package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.T;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Assert;
import org.junit.Test;

public class AmberSitesFileTest
{
    private RefGenomeVersion mGenomeVersion;
    ListMultimap<Chromosome, AmberSite> originalData;

    @Test
    public void writeReadTest() throws Exception
    {
        mGenomeVersion = RefGenomeVersion.V37;
        testForGenomeVersion();
        mGenomeVersion = RefGenomeVersion.V38;
        testForGenomeVersion();
    }

    private void testForGenomeVersion() throws IOException
    {
        originalData = ArrayListMultimap.create();
        putData(_1, 1000, A, C, true, 0.11);
        putData(_1, 1234, A, G, false, 0.56);
        putData(_1, 4564, T, G, false, 0.86453211);

        putData(_2, 232, T, A, false, 0.342);
        putData(_2, 232, T, G, true, 0.53211);

        putData(_X, 23782, T, G, true, 0.53211);
        putData(_X, 77999, T, A, false, 0.342);

        putData(_Y, 23782, T, G, true, 0.53211);
        putData(_X, 77999, T, A, false, 0.342);

        File tempDir = Files.createTempDirectory("amber").toFile();
        File outputFile = new File(tempDir, "ambersites.tsv.gz");
        AmberSitesFile.writeData(originalData, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, AmberSite> readData = AmberSitesFile.loadFile(outputFile.getAbsolutePath());

        // The original and new should be equal because equality for AmberSite
        // uses just position and bases. However, we need to check the snp check
        // flag and the frequency, which will have been truncated.
        Assert.assertEquals(originalData, readData);
        List<AmberSite> readData1 = readData.get(_1);
        Assert.assertEquals(originalData.get(_1).get(0).VariantAlleleFrequency, readData1.get(0).VariantAlleleFrequency, 0.0001);
        Assert.assertEquals(originalData.get(_1).get(0).snpCheck(), readData1.get(0).snpCheck());
        Assert.assertEquals(originalData.get(_1).get(2).VariantAlleleFrequency, readData1.get(2).VariantAlleleFrequency, 0.0001);
        Assert.assertEquals(originalData.get(_1).get(2).snpCheck(), readData1.get(2).snpCheck());
    }

    private void putData(HumanChromosome chromosome, int position, AmberBase ref, AmberBase alt, boolean snpCheck, double frequency)
    {
        originalData.put(chromosome, as(chromosome, position, ref, alt, snpCheck, frequency));
    }

    private AmberSite as(HumanChromosome chromosome, int position, AmberBase ref, AmberBase alt, boolean snpCheck, double frequency)
    {
        return new AmberSite(mGenomeVersion.versionedChromosome(chromosome), position, ref.name(), alt.name(), snpCheck, frequency);
    }
}
