package com.hartwig.hmftools.amber.e2e;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AmberScenario
{
    public static void main(String[] args) throws Exception
    {
        File refGenomeFile = new File("/Users/timlavers/work/data/reference_genomes/38/Homo_sapiens_assembly38.alt.masked.fasta");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
        File junkDir = new File("/Users/timlavers/work/junk");
        File outputFile = new File(junkDir, "NoBafs.tumor.bam");
        new AmberScenario("NoBafs").TumorReads.writeBam(refGenomeSource, outputFile.getAbsolutePath());
    }

    private final String Name;
    private final AmberReadsSpecification TumorReads;
    private final AmberReadsSpecification ReferenceReads;

    public AmberScenario(final String name) throws Exception
    {
        Name = name;
        var amberSites = AmberSitesFile.loadFile(getSitesFile().getAbsolutePath());
        File tumorReadsFile = getTumorReadsFile();
        if(tumorReadsFile.exists())
        {
            TumorReads = new AmberReadsSpecification(amberSites, tumorReadsFile);
        }
        else
        {
            TumorReads = null;
        }
        File referenceReadsFile = getReferenceReadsFile();
        if(referenceReadsFile.exists())
        {
            ReferenceReads = new AmberReadsSpecification(amberSites, referenceReadsFile);
        }
        else
        {
            ReferenceReads = null;
        }
        Preconditions.checkState(TumorReads != null || ReferenceReads != null);
    }

    public void checkResults(Multimap<Chromosome, AmberBAF> results)
    {
        Assert.assertEquals(TumorReads.chromosomes(), results.keySet());
        TumorReads.chromosomes().forEach(c -> checkTumorChromosomeResults(results.get(c), TumorReads.specifications(c)));
    }

    File getTumorBamFile()
    {
        if(TumorReads == null)
        {
            return null;
        }
        return new File(testResourcesDir(), Name + ".tumor.bam");
    }

    String getTumorSampleName()
    {
        return Name + "_tumor";
    }

    File createAmberLocationsFile(File destination) throws IOException
    {
        File sitesFile = getSitesFile();
        File locationsFile = new File(destination, "AmberGermlineSites.38.tsv");
        FileUtils.copyFile(sitesFile, locationsFile);
        return locationsFile;
    }

    private void checkTumorChromosomeResults(Collection<AmberBAF> results, List<AmberSiteExpectation> expectedResults)
    {
        Assert.assertEquals(results.size(), expectedResults.size());
        for(AmberBAF baf : results)
        {
            AmberSiteExpectation match = expectedResults.stream().filter(expectation -> expectation.aligns(baf)).findFirst().orElseThrow();
            match.checkTumorBaf(baf);
        }
    }

    private File testResourcesDir()
    {
        return Path.of("src", "test", "resources", "e2e").toFile();
    }

    private File getSitesFile()
    {
        return new File(testResourcesDir(), Name + ".tsv");
    }

    private File getTumorReadsFile()
    {
        return new File(testResourcesDir(), Name + ".tumor.reads.tsv");
    }

    private File getReferenceReadsFile()
    {
        return new File(testResourcesDir(), Name + ".reference.reads.tsv");
    }
}
