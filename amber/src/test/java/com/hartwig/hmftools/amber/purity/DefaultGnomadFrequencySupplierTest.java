package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.*;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Assert;
import org.junit.Test;

public class DefaultGnomadFrequencySupplierTest
{
    ListMultimap<Chromosome, AmberSite> sites = ArrayListMultimap.create();

    @Test
    public void getFrequencyTest()
    {
        addSite(_1, 1_000, 0.01);
        addSite(_1, 2_000, 0.02);
        addSite(_1, 3_000, 0.03);
        addSite(_1, 4_000, 0.04);
        addSite(_2, 2_000, 0.02);
        addSite(_3, 3_000, 0.24);
        DefaultGnomadFrequencySupplier supplier = new DefaultGnomadFrequencySupplier(sites, V38);
        Assert.assertEquals(0.01, supplier.getFrequency(V38.versionedChromosome(_1), 1000), 0.00001);
        Assert.assertEquals(0.02, supplier.getFrequency(V38.versionedChromosome(_1), 2000), 0.00001);
        Assert.assertEquals(0.24, supplier.getFrequency(V38.versionedChromosome(_3), 3000), 0.00001);
    }

    @Test
    public void throwExceptionForUnknownSiteTest()
    {
        DefaultGnomadFrequencySupplier supplier = new DefaultGnomadFrequencySupplier(sites, V38);
        Assert.assertThrows(IllegalArgumentException.class, () -> supplier.getFrequency(V38.versionedChromosome(_1), 1000));
    }

    private void addSite(HumanChromosome chromosome, int position, double gnomadFrequency)
    {
        sites.put(chromosome, as(chromosome, position, gnomadFrequency));
    }

    private AmberSite as(HumanChromosome chromosome, int position, double gnomadFrequency)
    {
        return new AmberSite(V38.versionedChromosome(chromosome), position, "A", "G", false, gnomadFrequency);
    }
}
