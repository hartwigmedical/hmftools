package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.*;

import static com.hartwig.hmftools.bamtools.remapper.RemapperTestBase.bwa;

public class AltContigRemapperTest extends RemapperTestBase
{

    @Test
    public void runTest()
    {
        Assert.assertEquals(14, records.size()); // Sanity check

        // Create a remapper that will remap the tiny bam test data file
        // using the HlaAligner (see below) test aligner (rather than one that uses the entire genome).
        File tempDir = new File("/Users/timlavers/work/junk"); // TODO is there a preferred way of handling temp files?
        File outputFile = new File(tempDir, "test.bam");
        File inputFile = getTestFile("tiny.bam");
        Aligner aligner = new HlaAligner(records);
        AltContigRemapperConfig config = new DummyConfig(aligner, inputFile.getAbsolutePath(), outputFile.getAbsolutePath());
        AltContigRemapper remapper = new AltContigRemapper(config);

        // Run the alignment.
        remapper.run();

        // Check that there are no HLA references in the dictionary of the output file;
        List<SAMSequenceRecord> dictionaryRecords = readDictionarySequences(outputFile);
        dictionaryRecords.forEach(record -> Assert.assertFalse(record.getSequenceName().contains("HLA")));

        // Read in the input file as a list of SAM records.
        List<SAMRecord> results = readSamFile(outputFile);

        // Two of the alignment sequences returned have 2 elements,
        // but two of the input sequences are supplementary HLA alignments so are ignored.
        Assert.assertEquals(records.size(), results.size());

        // The first four records in our test data file should not get altered and are already in order.
        check(records.get(0), results.get(0));
        check(records.get(1), results.get(1));
        check(records.get(2), results.get(2));
        check(records.get(3), results.get(3));

        // Expected[4] is input[4] with main and mate alignments remapped.
        // Main alignment: 	chr1, 194358863 -> chr1, 194358862 (see A00624:8:HHKYHDSXX:2:1516:12156:4225)
        // Mate alignment  HLA-A*01:03, 3188 -> chr6, 29945439.
        SAMRecord expected4 = records.get(4).deepCopy();
        expected4.setReferenceIndex(0);
        expected4.setAlignmentStart(194358862);
        expected4.setMappingQuality(60);
        expected4.setCigarString("35S116M");
        expected4.setMateAlignmentStart(29945439);
        expected4.setMateReferenceIndex(5);
        check(expected4, results.get(4));

        // Expected[5] is input[8] remapped
        // Main alignment: HLA-B*27:05:18, 3006 -> chr2, 32916241 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*40:06:01:01, 2533 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected5 = records.get(8).deepCopy();
        expected5.setReferenceIndex(1);
        expected5.setAlignmentStart(32916241);
        expected5.setMappingQuality(16);
        expected5.setCigarString("42S32M77S");
        expected5.setMateReferenceIndex(5);
        expected5.setMateAlignmentStart(31354760);
        check(expected5, results.get(5));

        // Expected[6] is input[13] remapped
        // Both the main alignment and the mate refer to HLA-B*08:19N
        // which gets remapped to chr2 for the main alignment and to
        // chr6 for the mate.
        // See the mappings for A00624:8:HHKYHDSXX:4:1304:7075:14998 below.
        SAMRecord expected6 = records.get(13).deepCopy();
        expected6.setReferenceIndex(1);
        expected6.setAlignmentStart(32916297);
        expected6.setMappingQuality(0);
        expected6.setCigarString("67S36M48S");
        expected6.setMateReferenceIndex(5);
        expected6.setMateAlignmentStart(31354513);
        check(expected6, results.get(6));

        // Expected[7] is the result of remapping input[7].
        // Main alignment: HLA-A*01:03, 3188 -> chr6, 29945439
        // Mate alignment: HLA-A*01:03, 3188 -> chr1, 194358862  (see A00624:8:HHKYHDSXX:2:1516:12156:4225)
        SAMRecord expected7 = records.get(7).deepCopy();
        expected7.setReferenceIndex(5);
        expected7.setAlignmentStart(29945439);
        expected7.setMappingQuality(60);
        expected7.setCigarString("151M");
        expected7.setMateReferenceIndex(0);
        expected7.setMateAlignmentStart(194358862);
        check(expected7, results.get(7));

        // Expected[8] is a second remapping of input[8]
        // Main alignment: HLA-B*27:05:18 -> chr6, 31354347 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*40:06:01:01 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected8 = records.get(8).deepCopy();
        expected8.setReferenceIndex(5);
        expected8.setAlignmentStart(31354347);
        expected8.setMappingQuality(60);
        expected8.setCigarString("72M79S");
        expected8.setMateReferenceIndex(5);
        expected8.setMateAlignmentStart(31354760);
        check(expected8, results.get(8));

        // Expected[9] is a second remapping of input[13]
        // Main alignment: HLA-B*08:19N -> chr6, 31354375 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        // Mate alignment: HLA-B*08:19N -> chr6, 31354513 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        SAMRecord expected9 = records.get(13).deepCopy();
        expected9.setReferenceIndex(5);
        expected9.setAlignmentStart(31354375);
        expected9.setMappingQuality(60);
        expected9.setCigarString("44M107S");
        expected9.setMateReferenceIndex(5);
        expected9.setMateAlignmentStart(31354513);
        check(expected9, results.get(9));

        // Element 10 is the result of remapping record 12.
        // Main alignment: HLA-B*08:19N, 2778 -> chr6, 31354513
        // Mate alignment: HLA-B*08:19N, 3023-> chr6, 31354375 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected10 = records.get(12).deepCopy();
        expected10.setReferenceIndex(5);
        expected10.setAlignmentStart(31354513);
        expected10.setMappingQuality(60);
        expected10.setCigarString("151M");
        expected10.setMateReferenceIndex(5);
        expected10.setMateAlignmentStart(31354375);
        check(expected10, results.get(10));

        // Element 11 is the result of remapping record 9.
        // Main alignment: HLA-B*40:06:01:01 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*27:05:18 -> chr6, 31354347 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected11 = records.get(9).deepCopy();
        expected11.setReferenceIndex(5);
        expected11.setAlignmentStart(31354760);
        expected11.setMappingQuality(60);
        expected11.setCigarString("151M");
        expected11.setMateReferenceIndex(5);
        expected11.setMateAlignmentStart(31354347);
        check(expected11, results.get(11));

        // Element 12 is the result of remapping record 11.
        // Main alignment: HLA-B*07:05:01 -> chr6, 31355729 (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        // Mate alignment: HLA-B*07:05:01 -> chr6, 31356297 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected12 = records.get(11).deepCopy();
        expected12.setReferenceIndex(5);
        expected12.setAlignmentStart(31355729);
        expected12.setMappingQuality(60);
        expected12.setCigarString("151M");
        expected12.setMateReferenceIndex(5);
        expected12.setMateAlignmentStart(31356297);
        check(expected12, results.get(12));

        // Element 13 is the result of remapping record 10.
        // Main alignment: HLA-B*07:05:01, 711 -> chr6, 31356297  (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        // Mate alignment: HLA-B*07:05:01, 1279 -> chr6, 31355729  (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        SAMRecord expected13 = records.get(10).deepCopy();
        expected13.setReferenceIndex(5);
        expected13.setAlignmentStart(31356297);
        expected13.setMappingQuality(60);
        expected13.setCigarString("151M");
        expected13.setMateReferenceIndex(5);
        expected13.setMateAlignmentStart(31355729);
        check(expected13, results.get(13));
    }
}

class DummyConfig extends AltContigRemapperConfig
{
    private final Aligner aligner;

    DummyConfig(Aligner aligner, String inputBam, String outputBam)
    {
        super(inputBam, outputBam, "/Users/timlavers/work/apps/samtools-1.21/samtools"); // TODO how should this be specified?
        this.aligner = aligner;
    }

    @Override
    public Aligner aligner()
    {
        return aligner;
    }
}

/**
 * This class implements Aligner by hard-coding the alignments for the hla alt records
 * in the test data file tiny.sam to the values returned by a BwaAligner instance that uses
 * the hg38_no_alts reference sequence.
 */
class HlaAligner implements Aligner
{

    private final List<SAMRecord> records;
    private final Map<Integer, List<BwaMemAlignment>> recordIndexToAlignments = new HashMap<>();

    HlaAligner(List<SAMRecord> records)
    {
        this.records = records;

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(4).getReadName()); // Sanity check
        recordIndexToAlignments.put(4, List.of(bwa("0,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,-1,-1,0")));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(7).getReadName()); // Sanity check
        recordIndexToAlignments.put(7, List.of(bwa("0,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,-1,-1,0")));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(8).getReadName()); // Sanity check
        recordIndexToAlignments.put(8, List.of(
                bwa("16,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,-1,-1,0"),
                bwa("2048,1,32916241,32916273,42,74,16,0,32,28,42S32M77S,32,xaTagIgnored,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(9).getReadName()); // Sanity check
        recordIndexToAlignments.put(9, List.of(
                bwa("16,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(10).getReadName()); // Sanity check
        recordIndexToAlignments.put(10, List.of(
                bwa("16,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(11).getReadName()); // Sanity check
        recordIndexToAlignments.put(11, List.of(
                bwa("16,5,31355729,31355880,0,151,60,0,151,57,151M,151,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(12).getReadName()); // Sanity check
        recordIndexToAlignments.put(12, List.of(
                bwa("16,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(13).getReadName()); // Sanity check
        recordIndexToAlignments.put(13, List.of(
                bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0"),
                bwa("2048,1,32916297,32916333,67,103,0,0,36,36,67S36M48S,36,null,-1,-1,0")
        ));
    }

    @Override
    public List<BwaMemAlignment> alignSequence(byte[] bases)
    {
        for(int index = 0; index < records.size(); index++)
        {
            if(Arrays.equals(bases, records.get(index).getReadBases()))
            {
                return recordIndexToAlignments.get(index);
            }
        }
        throw new IllegalArgumentException("No alignment found for " + Arrays.toString(bases));
    }

    @Override
    public ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2)
    {
        List<BwaMemAlignment> resultLeft = alignSequence(bases1);
        List<BwaMemAlignment> resultRight = alignSequence(bases2);
        return new ImmutablePair<>(resultLeft, resultRight);
    }
}
