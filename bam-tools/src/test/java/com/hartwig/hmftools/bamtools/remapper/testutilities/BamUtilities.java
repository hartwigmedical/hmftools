package com.hartwig.hmftools.bamtools.remapper.testutilities;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamUtilities
{
    public static void main(String[] args) throws Exception
    {
        produceBam();
    }

    public static void produceReducedChrFile() throws IOException
    {
        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        RefGenomeSource refGenomeSource =
                new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));
        int start = 0; //
        int end = start + 1_000_000;
        var chr = refGenomeSource.getBaseString("chr1", start, end);
        //        System.out.println(chr.substring(10000, 10100));
        File chrFile = new File(outputDir, "chr1_part_0.txt");
        Files.writeString(chrFile.toPath(), chr, StandardCharsets.UTF_8);
    }

    public static void produceBam() throws Exception
    {
//        loadAlignerLibrary(null);

        File outputDir = new File("/Users/timlavers/work/junk/rubbish");
        RefGenomeSource refGenomeSource =
                new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));
//        int start = 0; // 0 is \n

        var outputFileName = new File(outputDir, "bu5.bam").getAbsolutePath();
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

        RegionDepth chr1Depth = new RegionDepth(0, 100_000, 200_000, 100);
        chr1Depth.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource);
        bamWriter.close();
    }
}

