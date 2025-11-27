package com.hartwig.hmftools.amber.e2e;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.apache.commons.io.FileUtils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AmberScenario
{
    public static void main(String[] args) throws Exception
    {
        File refGenomeFile = new File("/Users/timlavers/work/data/reference_genomes/38/Homo_sapiens_assembly38.alt.masked.fasta");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
        File junkDir = new File("/Users/timlavers/work/junk");
        File outputFile = new File(junkDir, "TwoChromosomes.tumor.bam");
        new AmberScenario("TwoChromosomes").TumorReads.writeBam(refGenomeSource, outputFile.getAbsolutePath());
    }

    private final String Name;
    private final AmberReadsSpecification TumorReads;
    private final AmberReadsSpecification ReferenceReads;
    private final ListMultimap<Chromosome, AmberSiteRead> Sites = ArrayListMultimap.create();

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

    public File getTumorBamFile()
    {
        if (TumorReads == null)
        {
            return null;
        }
        return new File(testResourcesDir(), Name + ".tumor.bam");
    }

    public File createAmberLocationsFile(File destination) throws IOException
    {
        File sitesFile = getSitesFile();
        File locationsFile = new File(destination, "AmberGermlineSites.38.tsv");
        FileUtils.copyFile(sitesFile, locationsFile);
        return locationsFile;
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
