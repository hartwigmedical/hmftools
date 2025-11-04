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
//        File refGenomeFile =
//                new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
//        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));

       BamRecipe bamRecipe = new BamRecipe(new ConstantChromosomeLengths(5_000));
        int regionOffset = 1;
        ChromosomeRegionDepths depths = new GCRatioChromosomeRegionDepths(_1, 40);
        depths.addRange(regionOffset, regionOffset + 800, 1);
        bamRecipe.add(depths);
        regionOffset += 1000;

        depths = new GCRatioChromosomeRegionDepths(_1, 42);
        depths.addRange(regionOffset, regionOffset + 800, 1);
        bamRecipe.add(depths);
        regionOffset += 1000;

        depths = new GCRatioChromosomeRegionDepths(_1, 44);
        depths.addRange(regionOffset, regionOffset + 1000, 10);
        bamRecipe.add(depths);
        regionOffset += 1000;

        depths = new GCRatioChromosomeRegionDepths(_1, 46);
        depths.addRange(regionOffset, regionOffset + 1000, 10);
        bamRecipe.add(depths);
        regionOffset += 1000;

        depths = new GCRatioChromosomeRegionDepths(_1, 48);
        depths.addRange(regionOffset, regionOffset + 1000, 10);
        bamRecipe.add(depths);

        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(outputDir, "low_depth_regions.bam");
        bamRecipe.writeToBam(bamFile.getAbsolutePath());
    }
}

