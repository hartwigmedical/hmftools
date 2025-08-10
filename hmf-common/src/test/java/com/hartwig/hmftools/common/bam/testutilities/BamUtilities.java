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
        BamRecipe bamRecipe = new BamRecipe(new ConstantChromosomeLengths(101_000));
        int regionOffset = 1_001;
        //        BamRecipe bamRecipe = new BamRecipe(new RefGenomeBackedChromosomeLengths(refGenomeSource));
            for (int i=0; i<11; i++)
            {
                double ratio = 1.0 * (40 + i) / 100.0;
                ChromosomeRegionDepths depths = new GCRatioChromosomeRegionDepths(0, ratio);
                depths.addRange(regionOffset, regionOffset + 1_000, 10 +i);
                regionOffset += 1000;
                depths.addRange(regionOffset, regionOffset + 1_000, 10 +i);
                regionOffset += 1000;
                depths.addRange(regionOffset, regionOffset + 1_000, 10 +i);
                regionOffset += 1000;
                depths.addRange(regionOffset, regionOffset + 1_000, 20 +i);
                regionOffset += 1000;
                depths.addRange(regionOffset, regionOffset + 1_000, 50 +i);
                regionOffset += 1000;
                bamRecipe.add(depths);
            }
//        for (int i=0; i<100; i++)
//        {
//            double ratio = 1.0 * i / 100.0;
//            ChromosomeRegionDepths depths = new GCRatioChromosomeRegionDepths(0, ratio);
//            depths.addRange(regionOffset, regionOffset + 1_000, 1);
//            regionOffset += 1000;
//            bamRecipe.add(depths);
//        }
//        ChromosomeRegionDepths chr1Depths = new GCRatioChromosomeRegionDepths(0, 0.20);
//        chr1Depths.addRange(regionOffset, regionOffset + 1_000, 10);
//        chr1Depths.addRange(regionOffset + 1000, regionOffset + 2_000, 100);
//        bamRecipe.add(chr1Depths);
//        chr1Depths.addRange(regionOffset + 2000, regionOffset + 3_000, 50);
        //        chr1Depths.addRange(regionOffset + 3_000, regionOffset + 6_000, 10);
        //        chr1Depths.addRange(regionOffset + 6_000, regionOffset + 9_000, 100);

//        ChromosomeRegionDepths chrXDepths = new RepeatingACGTChromosomeRegionDepths(22);
//        regionOffset = 1_001;
//        chrXDepths.addRange(regionOffset, regionOffset + 3_000, 10);
//        chr1Depths.addRange(regionOffset + 3_000, regionOffset + 6_000, 10);
//        chr1Depths.addRange(regionOffset + 6_000, regionOffset + 9_000, 100);
//        bamRecipe.add(chrXDepths);

//        ChromosomeRegionDepths chr2Depths = new GenomeBackedChromosomeRegionDepths(1, refGenomeSource);
//        chr2Depths.addRange(regionOffset, regionOffset + 100_000, 50);
//        chr2Depths.addRange(regionOffset + 200_000, regionOffset + 300_000, 50);
//        chr2Depths.addRange(regionOffset + 400_000, regionOffset + 500_000, 50);
//        bamRecipe.add(chr2Depths);

        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(outputDir, "Example10.bam");
        bamRecipe.writeToBam(bamFile.getAbsolutePath());
    }
}

