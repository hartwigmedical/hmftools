package com.hartwig.hmftools.common.bam.testutilities;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.*;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import java.io.File;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.apache.commons.lang3.RandomUtils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamMaker
{
    public static void main(String[] args) throws Exception
    {
        File refGenomeFile =
                new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));

       BamRecipe bamRecipe = new BamRecipe(new ConstantChromosomeLengths(101_000));
        int regionOffset = 1;
//        BamRecipe bamRecipe = new BamRecipe(new RefGenomeBackedChromosomeLengths(refGenomeSource));
//        ChromosomeRegionDepths chr1Depths = new GCRatioChromosomeRegionDepths(_1, 0.5);
//        ChromosomeRegionDepths chr2Depths = new GCRatioChromosomeRegionDepths(_2, 0.5);
//        ChromosomeRegionDepths chrXDepths = new GCRatioChromosomeRegionDepths(_X, 0.5);
//        ChromosomeRegionDepths chrYDepths = new GCRatioChromosomeRegionDepths(_Y, 0.5);
//        chr1Depths.addRange(regionOffset, regionOffset + 1000, 100);
//        chr2Depths.addRange(regionOffset, regionOffset + 1000, 110);
//        chrXDepths.addRange(regionOffset, regionOffset + 1000, 60);
//        chrYDepths.addRange(regionOffset, regionOffset + 1000, 40);
//        regionOffset += 1000;
//        chr1Depths.addRange(regionOffset, regionOffset + 1000, 110);
//        chr2Depths.addRange(regionOffset, regionOffset + 1000, 100);
//        chrXDepths.addRange(regionOffset, regionOffset + 1000, 60);
//        chrYDepths.addRange(regionOffset, regionOffset + 1000, 40);
//        regionOffset += 1000;
//        chr1Depths.addRange(regionOffset, regionOffset + 1000, 104);
//        chr2Depths.addRange(regionOffset, regionOffset + 1000, 104);
//        chrXDepths.addRange(regionOffset, regionOffset + 1000, 50);
//        chrYDepths.addRange(regionOffset, regionOffset + 1000, 50);

//        bamRecipe.add(chr1Depths);
//        bamRecipe.add(chr2Depths);
//        bamRecipe.add(chrXDepths);
//        bamRecipe.add(chrYDepths);

        /*
        ChromosomeRegionDepths chr1Depths = new GCRatioChromosomeRegionDepths(_1, 0.5);
        ChromosomeRegionDepths chr2Depths = new GCRatioChromosomeRegionDepths(_2, 0.5);
        for(int i = 0; i <= 40; i++)
        {
            int random = RandomUtils.nextInt(0, 10);
            chr1Depths.addRange(regionOffset, regionOffset + 1_000, 95 + random);
            random = RandomUtils.nextInt(0, 100);
            chr2Depths.addRange(regionOffset, regionOffset + 1_000, 950 + random);
            regionOffset += 1000;
        }
        for(int i = 41; i <= 70; i++)
        {
            int random = RandomUtils.nextInt(0, 100);
            chr1Depths.addRange(regionOffset, regionOffset + 1_000, 950 + random);
            random = RandomUtils.nextInt(0, 10);
            chr2Depths.addRange(regionOffset, regionOffset + 1_000, 95 + random);
            regionOffset += 1000;
        }
        for(int i = 71; i <= 101; i++)
        {
            int random = RandomUtils.nextInt(0, 10);
            chr1Depths.addRange(regionOffset, regionOffset + 1_000, 95 + random);
            random = RandomUtils.nextInt(0, 100);
            chr2Depths.addRange(regionOffset, regionOffset + 1_000, 950 + random);
            regionOffset += 1000;
        }
        bamRecipe.add(chr1Depths);
        bamRecipe.add(chr2Depths);
        */
                for (int i=0; i<=100; i++)
                {
                    ChromosomeRegionDepths depths = new GCRatioChromosomeRegionDepths(_1, i);
                    depths.addRange(regionOffset, regionOffset + 1_000, 10);
                    regionOffset += 1000;
                    bamRecipe.add(depths);
                }
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

//                ChromosomeRegionDepths chr2Depths = new GenomeBackedChromosomeRegionDepths(_2, refGenomeSource);
//                chr2Depths.addRange(regionOffset, regionOffset + 100_000, 50);
//                chr2Depths.addRange(regionOffset + 200_000, regionOffset + 300_000, 50);
//                chr2Depths.addRange(regionOffset + 400_000, regionOffset + 500_000, 50);
//                bamRecipe.add(chr2Depths);

        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(outputDir, "Example14.bam");
        bamRecipe.writeToBam(bamFile.getAbsolutePath());
    }
}

