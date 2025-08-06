package com.hartwig.hmftools.common.bam.testutilities;

import java.io.File;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamUtilities
{
    public static void main(String[] args) throws Exception
    {
        File refGenomeFile = new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
        BamRecipe bamRecipe = new BamRecipe(new ConstantChromosomeLengths(5_000));
//        BamRecipe bamRecipe = new BamRecipe(new RefGenomeBackedChromosomeLengths(refGenomeSource));
        ChromosomeRegionDepths chr1Depths = new RepeatingACGTChromosomeRegionDepths(0);
        int regionOffset = 1_001;
        chr1Depths.addRange(regionOffset, regionOffset + 3_000, 100);
//        chr1Depths.addRange(regionOffset + 3_000, regionOffset + 6_000, 10);
//        chr1Depths.addRange(regionOffset + 6_000, regionOffset + 9_000, 100);
        bamRecipe.add(chr1Depths);

//        ChromosomeRegionDepths chr2Depths = new GenomeBackedChromosomeRegionDepths(1, refGenomeSource);
//        chr2Depths.addRange(regionOffset, regionOffset + 100_000, 50);
//        chr2Depths.addRange(regionOffset + 200_000, regionOffset + 300_000, 50);
//        chr2Depths.addRange(regionOffset + 400_000, regionOffset + 500_000, 50);
//        bamRecipe.add(chr2Depths);

        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(outputDir, "Example9.bam");
        bamRecipe.writeToBam(bamFile.getAbsolutePath());
    }
}

