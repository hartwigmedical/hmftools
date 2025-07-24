package com.hartwig.hmftools.common.bam.testutilities;

import java.io.File;

import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamUtilities
{
    public static void main(String[] args) throws Exception
    {
        BamRecipe bamRecipe = new BamRecipe();
        ChromosomeRegionDepths chr1Depths = new ChromosomeRegionDepths(0);
        int regionOffset = 10_000_000;
        chr1Depths.addRange(regionOffset, regionOffset + 100_000, 100);
        chr1Depths.addRange(regionOffset + 200_000, regionOffset + 300_000, 20);
        chr1Depths.addRange(regionOffset + 400_000, regionOffset + 500_000, 100);
        bamRecipe.add(chr1Depths);

        ChromosomeRegionDepths chr2Depths = new ChromosomeRegionDepths(1);
        chr2Depths.addRange(regionOffset, regionOffset + 100_000, 100);
        chr2Depths.addRange(regionOffset + 200_000, regionOffset + 300_000, 20);
        chr2Depths.addRange(regionOffset + 400_000, regionOffset + 500_000, 100);
        bamRecipe.add(chr2Depths);

        File refGenomeFile = new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(outputDir, "Example1.bam");
        bamRecipe.writeToBam(bamFile.getAbsolutePath(), refGenomeSource);
    }

    public static void produceBam() throws Exception
    {
        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        RefGenomeSource refGenomeSource =
                new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));
//        int start = 0; // 0 is \n

        var outputFileName = new File(outputDir, "T2.bam").getAbsolutePath();
        var header = new SAMFileHeader();
        header.getSequenceDictionary().addSequence(new SAMSequenceRecord("chr1", 248956422));
        SAMFileWriter bamWriter = new FastBamWriter(header, outputFileName);

//        ChromosomeRegionDepths chr1Depths = new ChromosomeRegionDepths(0);
//        chr1Depths.addRange(100_000, 100_200, 10);
//        chr1Depths.addRange(100_200, 100_400, 5);
//        chr1Depths.addRange(100_400, 100_600, 10);
//        chr1Depths.addRange(110_000, 119_999, 5);
//        chr1Depths.addRange(120_000, 129_999, 10);
//        chr1Depths.writeToBam(bamWriter, refGenomeSource);

        int regionOffset = 10_000_000;
        RegionDepth chr1Depth0 = new RegionDepth(0, regionOffset + 10_000, regionOffset + 20_000, 1000);
        chr1Depth0.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        RegionDepth chr1Depth1 = new RegionDepth(0, regionOffset + 20_000, regionOffset + 30_000, 100);
        chr1Depth1.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        RegionDepth chr1Depth2 = new RegionDepth(0, regionOffset + 30_000, regionOffset + 40_000, 1000);
        chr1Depth2.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        bamWriter.close();
    }

    public static void writeBam(BamRecipe bamRecipe, File bamFile, File refGenomeFile) throws Exception
    {
        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(refGenomeFile));
//        refGenomeSource.

        var outputFileName = bamFile.getAbsolutePath();
        var header = new SAMFileHeader();
        header.getSequenceDictionary().addSequence(new SAMSequenceRecord("chr1", 248956422));
        SAMFileWriter bamWriter = new FastBamWriter(header, outputFileName);

//        ChromosomeRegionDepths chr1Depths = new ChromosomeRegionDepths(0);
//        chr1Depths.addRange(100_000, 100_200, 10);
//        chr1Depths.addRange(100_200, 100_400, 5);
//        chr1Depths.addRange(100_400, 100_600, 10);
//        chr1Depths.addRange(110_000, 119_999, 5);
//        chr1Depths.addRange(120_000, 129_999, 10);
//        chr1Depths.writeToBam(bamWriter, refGenomeSource);

        int regionOffset = 10_000_000;
        RegionDepth chr1Depth0 = new RegionDepth(0, regionOffset + 10_000, regionOffset + 20_000, 1000);
        chr1Depth0.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        RegionDepth chr1Depth1 = new RegionDepth(0, regionOffset + 20_000, regionOffset + 30_000, 100);
        chr1Depth1.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        RegionDepth chr1Depth2 = new RegionDepth(0, regionOffset + 30_000, regionOffset + 40_000, 1000);
        chr1Depth2.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        bamWriter.close();
    }
}

