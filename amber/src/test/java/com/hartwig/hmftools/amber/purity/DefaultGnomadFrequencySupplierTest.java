package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.*;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

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
    public void throwExceptionForUnknownPositionTest()
    {
        DefaultGnomadFrequencySupplier supplier = new DefaultGnomadFrequencySupplier(sites, V38);
        Assert.assertThrows(IllegalArgumentException.class, () -> supplier.getFrequency(V38.versionedChromosome(_1), 1000));
    }

    @Test
    public void throwExceptionForUnknownChromosomeTest()
    {
        DefaultGnomadFrequencySupplier supplier = new DefaultGnomadFrequencySupplier(sites, V38);
        Assert.assertThrows(IllegalArgumentException.class, () -> supplier.getFrequency(V38.versionedChromosome(_Y), 1000));
    }

    //    @Test
    public void loadTest() throws Exception
    {
        sites = AmberSitesFile.loadFile("/Users/timlavers/work/junk/AmberGermlineSites.new.V38.tsv.gz");
        long before = getMemoryUsageEstimate();

        Stopwatch creationTimer = Stopwatch.createStarted();
        GnomadFrequencySupplier supplier = new DefaultGnomadFrequencySupplier(sites, V38);
        creationTimer.stop();
        System.out.println("creation time: " + creationTimer.elapsed(TimeUnit.MILLISECONDS));

        long after = getMemoryUsageEstimate();
        System.out.println("Memory: " + (after - before) / (1024 * 1024) + " MB");

        List<GenomePosition> sitesList = new ArrayList<>(sites.values());
        Stopwatch timer = Stopwatch.createStarted();
        for(GenomePosition site : sitesList)
        {
            supplier.getFrequency(site.chromosome(), site.position());
        }
        timer.stop();
        System.out.println("time: " + timer.elapsed(TimeUnit.MILLISECONDS));
    }

    private static long getMemoryUsageEstimate() throws InterruptedException
    {
        System.gc();
        Thread.sleep(100);
        System.gc();
        Thread.sleep(100);
        return Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
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
