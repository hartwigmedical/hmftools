package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.codon.Nucleotides;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.*;

import static com.hartwig.hmftools.bamtools.remapper.RemapperTestBase.bwa;
import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;

public class AltContigRemapperTest extends RemapperTestBase
{

    @Test
    public void runTest()
    {
        Assert.assertEquals(18, records.size()); // Sanity check

        // Create a remapper that will remap the tiny bam test data file
        // using the HlaAligner (see below) test aligner (rather than one that uses the entire genome).
        File tempDir = FileUtils.getTempDirectory();
        File outputFile = new File(tempDir, "test.bam");
        File inputFile = getTestFile("tiny.bam");
        PairAligner aligner = new HlaAligner(records);
        AltContigRemapperConfig config = new DummyConfig(aligner, inputFile.getAbsolutePath(), outputFile.getAbsolutePath());
        AltContigRemapper remapper = new AltContigRemapper(config);

        // Run the alignment.
        remapper.run();

        // Check that there are no HLA references in the dictionary of the output file;
        List<SAMSequenceRecord> dictionaryRecords = readDictionarySequences(outputFile);
        dictionaryRecords.forEach(record -> Assert.assertFalse(record.getSequenceName().contains("HLA")));

        // Read in the input file as a list of SAM records.
        List<SAMRecord> results = readSamFile(outputFile);

        // Two of the alignment sequences returned have 2 elements.
        // For both of theses, the input record produces two output records.
        // Two of the input records are HLA but are supplementary,
        // so are ignored.
        // The nett result is that the number of outputs is the number of inputs.
        Assert.assertEquals(records.size(), results.size());

        // In this unit test we do not set the samtools path,
        // so the output file is not sorted. This means that
        // the order of outputs is the order of inputs, taking
        // into account those records that produce either zero
        // or two outputs, as mentioned above.
        // The first four records in our test data file should not get altered.
        checkOutputContains(records.get(0), results);
        checkOutputContains(records.get(1), results);
        checkOutputContains(records.get(2), results);
        checkOutputContains(records.get(3), results);

        // input[4]: A00624:8:HHKYHDSXX:2:1516:12156:4225, reverse strand
        // BWA returns 81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0
        SAMRecord expectedFor4 = records.get(4).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", expectedFor4.getReadName()); // sanity check
        Assert.assertEquals(81, expectedFor4.getFlags()); // sanity check
        expectedFor4.setFlags(81);
        expectedFor4.setReferenceIndex(0);
        expectedFor4.setAlignmentStart(194358862 + 1);
        expectedFor4.setMappingQuality(60);
        expectedFor4.setCigarString("35S116M");
        expectedFor4.setMateAlignmentStart(29945439 + 1);
        expectedFor4.setMateReferenceIndex(5);
        expectedFor4.setInferredInsertSize(0);
        checkOutputContains(expectedFor4, results);

        // input[5] and input[6] produce no output.

        // input[7]: A00624:8:HHKYHDSXX:2:1516:12156:4225, forward strand
        // BWA returns 161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0
        SAMRecord expectedFor7 = records.get(7).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", expectedFor7.getReadName()); // sanity check
        Assert.assertEquals(161, expectedFor7.getFlags()); // sanity check
        expectedFor7.setFlags(161);
        expectedFor7.setReferenceIndex(5);
        expectedFor7.setAlignmentStart(29945439 + 1);
        expectedFor7.setMappingQuality(60);
        expectedFor7.setCigarString("151M");
        expectedFor7.setMateReferenceIndex(0);
        expectedFor7.setMateAlignmentStart(194358862 + 1);
        expectedFor7.setInferredInsertSize(0);
        checkOutputContains(expectedFor7, results);

        // input[8]: A00624:8:HHKYHDSXX:1:2559:3224:21292, reverse strand
        // BWA returns 81,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,5,31354347,-564
        SAMRecord expectedMainFor8 = records.get(8).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", expectedMainFor8.getReadName()); // sanity check
        Assert.assertEquals(145, expectedMainFor8.getFlags()); // sanity check
        expectedMainFor8.setFlags(81 + 2); // bumped by proper pair flag (2)
        expectedMainFor8.setReferenceIndex(5);
        expectedMainFor8.setAlignmentStart(31354760 + 1);
        expectedMainFor8.setMappingQuality(60);
        expectedMainFor8.setCigarString("151M");
        expectedMainFor8.setMateReferenceIndex(5);
        expectedMainFor8.setMateAlignmentStart(31354347 + 1);
        expectedMainFor8.setInferredInsertSize(-564);
        checkOutputContains(expectedMainFor8, results);

        // input[8]: A00624:8:HHKYHDSXX:1:2559:3224:21292, forward strand
        // BWA returns:
        // principal alignment: 161,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,5,31354760,564
        // secondary alignment: 2225,1,32916241,32916273,42,74,13,0,32,28,42S32M77S,32,null,5,31354760,0
        // principal:
        SAMRecord expectedSupplementaryFor8 = records.get(9).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", expectedSupplementaryFor8.getReadName()); // sanity check
        Assert.assertEquals(97, expectedSupplementaryFor8.getFlags()); // sanity check
        expectedSupplementaryFor8.setFlags(161 + 2); // gets proper pair flag
        expectedSupplementaryFor8.setReferenceIndex(5);
        expectedSupplementaryFor8.setAlignmentStart(31354347 + 1);
        expectedSupplementaryFor8.setMappingQuality(60);
        expectedSupplementaryFor8.setCigarString("72M79S");
        expectedSupplementaryFor8.setMateReferenceIndex(5);
        expectedSupplementaryFor8.setMateAlignmentStart(31354760 + 1);
        expectedSupplementaryFor8.setInferredInsertSize(564);
        checkOutputContains(expectedSupplementaryFor8, results);

        // secondary:
        SAMRecord expectedFor9 = records.get(9).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", expectedFor9.getReadName()); // sanity check
        Assert.assertEquals(97, expectedFor9.getFlags()); // sanity check
        expectedFor9.setFlags(2225);
        expectedFor9.setReferenceIndex(1);
        expectedFor9.setAlignmentStart(32916241 + 1);
        expectedFor9.setMappingQuality(13);
        expectedFor9.setCigarString("42S32M77S");
        expectedFor9.setMateReferenceIndex(5);
        expectedFor9.setMateAlignmentStart(31354760 + 1);
        expectedFor9.setInferredInsertSize(0);
        // This is a reverse strand match.
        expectedFor9.setReadBases(Nucleotides.reverseComplementBases(expectedFor9.getReadBases()));
        expectedFor9.setBaseQualities(reverseArray(expectedFor9.getBaseQualities()));
        checkOutputContains(expectedFor9, results);

        // input[10]: A00624:8:HHKYHDSXX:1:1446:18213:29684, forward strand
        // BWA returns 97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719
        SAMRecord expectedFor10 = records.get(10).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", expectedFor10.getReadName()); // sanity check
        Assert.assertEquals(163, expectedFor10.getFlags()); // sanity check
        expectedFor10.setFlags(97 + 2); // proper pair bump
        expectedFor10.setReferenceIndex(5);
        expectedFor10.setAlignmentStart(31355729 + 1);
        expectedFor10.setMappingQuality(60);
        expectedFor10.setCigarString("151M");
        expectedFor10.setMateReferenceIndex(5);
        expectedFor10.setMateAlignmentStart(31356297 + 1);
        expectedFor10.setInferredInsertSize(719);
        checkOutputContains(expectedFor10, results);

        // input[11]: A00624:8:HHKYHDSXX:1:1446:18213:29684, reverse strand
        // BWA returns 145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719
        SAMRecord expectedFor11 = records.get(11).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", expectedFor11.getReadName()); // sanity check
        Assert.assertEquals(83, expectedFor11.getFlags()); // sanity check
        expectedFor11.setFlags(145 + 2); // proper pair
        expectedFor11.setReferenceIndex(5);
        expectedFor11.setAlignmentStart(31356297 + 1);
        expectedFor11.setMappingQuality(60);
        expectedFor11.setCigarString("151M");
        expectedFor11.setMateReferenceIndex(5);
        expectedFor11.setMateAlignmentStart(31355729 + 1);
        expectedFor11.setInferredInsertSize(-719);
        checkOutputContains(expectedFor11, results);

        /*
        The last two inputs result in three lines of output, but in a
        pretty confusing way. The input lines are:
        A00624:8:HHKYHDSXX:4:1304:7075:14998	163	HLA-B*08:19N	2778	0	151M	=	3023	289	TGGAAAAGGAGGGAGCTACTCTCAGGCTGCGTGTAAGTGATGGGGGTGGGAGTGTGGAGGAGCTCACCCACCCCATAATTCCTCCTGTCCCACGTCTCCTGCGGGCTCTGACCAGGTCCTGTTTTTGTTCTACTCCAGGCAGCGACAGTGC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:151	MC:Z:107S44M	MQ:i:0	AS:i:151	XS:i:151
        A00624:8:HHKYHDSXX:4:1304:7075:14998	83	HLA-B*08:19N	3023	0	107S44M	=	2778	-289	GCGGGCAGTGGCCCGGGGTTTTTTTTGTTCACAAACGGTTTTAATTGGGTGTGTTTGGGGGGTTTGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCTGGGGGGAAAGGCCTGGGTAATGGAGATTCTTTGATTGGGATGTTTC	,,,,,,,,,,,,:,,,F:,,,,:,,,,,,,,,,:,,,,,:F,,F,:,,,,,,,,,,,F:F:,,,,,,F:FFFFFFFFFFFFFFF,,,FFFFF:F:F:F:FFF,,,:,FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:44	MC:Z:151M	MQ:i:0	AS:i:44	XS:i:48

        In an actual direct mapping without alts, the corresponding reads get mapped to:
        A00624:8:HHKYHDSXX:4:1304:7075:14998	2163	chr2	32916352	0	67H36M48H	chr6	31354514	0	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	F:FFFFFFFFFFFFFFF,,,FFFFF:F:F:F:FFF,	NM:i:0	MD:Z:36	MC:Z:151M	MQ:i:60	AS:i:36	XS:i:36	SA:Z:chr6,31354376,+,44M107S,60,1;
        A00624:8:HHKYHDSXX:4:1304:7075:14998	99	chr6	31354376	60	44M107S	=	31354514	289	GAAACATCCCAATCAAAGAATCTCCATTACCCAGGCCTTTCCCCCCAGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACAAACCCCCCAAACACACCCAATTAAAACCGTTTGTGAACAAAAAAAACCCCGGGCCACTGCCCGC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF,:,,,FFF:F:F:F:FFFFF,,,FFFFFFFFFFFFFFF:F,,,,,,:F:F,,,,,,,,,,,:,F,,F:,,,,,:,,,,,,,,,,:,,,,:F,,,:,,,,,,,,,,,,	NM:i:1	MD:Z:22C21	MC:Z:151M	MQ:i:60	AS:i:39	XS:i:20	SA:Z:chr2,32916352,-,67S36M48S,0,0;
        A00624:8:HHKYHDSXX:4:1304:7075:14998	147	chr6	31354514	60	151M	=	31354376	-289	GCACTGTCGCTGCCTGGAGTAGAACAAAAACAGGACCTGGTCAGAGCCCGCAGGAGACGTGGGACAGGAGGAATTATGGGGTGGGTGAGCTCCTCCACACTCCCACCCCCATCACTTACACGCAGCCTGAGAGTAGCTCCCTCCTTTTCCA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	MD:Z:111C39	MC:Z:44M107S	MQ:i:60	AS:i:146	XS:i:84

        The first of these is a hard-clipped supplementary alignment corresponding to the second alt alignment read, reverse strand.
        The second is a forward-strand version of the second alt alignment read.
        The third is a reverse-strand version of the first alt alignment read.
         */
        // input[12]: A00624:8:HHKYHDSXX:4:1304:7075:14998, forward strand
        // BWA returns 145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289
        SAMRecord expectedFor12 = records.get(12).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", expectedFor12.getReadName()); // sanity check
        Assert.assertEquals(163, expectedFor12.getFlags()); // sanity check
        expectedFor12.setFlags(145 + 2); // proper pair
        expectedFor12.setReferenceIndex(5);
        expectedFor12.setAlignmentStart(31354513 + 1);
        expectedFor12.setMappingQuality(60);
        expectedFor12.setCigarString("151M");
        expectedFor12.setMateReferenceIndex(5);
        expectedFor12.setMateAlignmentStart(31354375 + 1);
        expectedFor12.setInferredInsertSize(-289);
        // This is a reverse strand match.
        expectedFor12.setReadBases(Nucleotides.reverseComplementBases(expectedFor12.getReadBases()));
        expectedFor12.setBaseQualities(reverseArray(expectedFor12.getBaseQualities()));
        checkOutputContains(expectedFor12, results);

        // input[13]: A00624:8:HHKYHDSXX:4:1304:7075:14998, reverse strand
        // BWA returns a supplementary forward strand alignment and a reverse strand principal alignment:
        // supplementary: 2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0
        // principal: 97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289

        SAMRecord expectedPrincipalFor13 = records.get(13).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", expectedPrincipalFor13.getReadName()); // sanity check
        Assert.assertEquals(83, expectedPrincipalFor13.getFlags()); // sanity check
        expectedPrincipalFor13.setFlags(97 + 2); // proper pair
        expectedPrincipalFor13.setReferenceIndex(5);
        expectedPrincipalFor13.setAlignmentStart(31354375 + 1);
        expectedPrincipalFor13.setMappingQuality(60);
        expectedPrincipalFor13.setCigarString("44M107S");
        expectedPrincipalFor13.setMateReferenceIndex(5);
        expectedPrincipalFor13.setMateAlignmentStart(31354513 + 1);
        expectedPrincipalFor13.setInferredInsertSize(289);
        // This is a forward strand match, whereas the input was reverse strand.
        expectedPrincipalFor13.setReadBases(Nucleotides.reverseComplementBases(expectedPrincipalFor13.getReadBases()));
        expectedPrincipalFor13.setBaseQualities(reverseArray(expectedPrincipalFor13.getBaseQualities()));
        checkOutputContains(expectedPrincipalFor13, results);

        SAMRecord expectedSecondaryFor13 = records.get(13).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", expectedSecondaryFor13.getReadName()); // sanity check
        Assert.assertEquals(83, expectedSecondaryFor13.getFlags()); // sanity check
        expectedSecondaryFor13.setFlags(2161);
        expectedSecondaryFor13.setReferenceIndex(1);
        expectedSecondaryFor13.setAlignmentStart(32916486 + 1);
        expectedSecondaryFor13.setMappingQuality(0);
        expectedSecondaryFor13.setCigarString("67S36M48S");
        expectedSecondaryFor13.setMateReferenceIndex(5);
        expectedSecondaryFor13.setMateAlignmentStart(31354513 + 1);
        expectedSecondaryFor13.setInferredInsertSize(0);
        checkOutputContains(expectedSecondaryFor13, results);

        // input[14] and input[15]: A00624:8:HHKYHDSXX:2:2276:27299:7623, forward and reverse respectively
        // These are mapped in the input to HLA-C. The data returned by BWA is:
        // 97,5,31356204,31356355,0,151,19,4,131,121,151M,21T0C19C79G28,"chr6,+31271111,151M,6;",5,31271513,-84542 (forward)
        // 145,5,31271513,31271664,0,151,18,6,121,111,151M,5G8C0T28G42T51C11,"chr6,-31356606,48M3I100M,10;",5,31356204,84542 (reverse)
        // These are not concordant (note the large insertion lengths) so an alt
        // alignment is used for the forward strand.
        SAMRecord expectedFor14 = records.get(14).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", expectedFor14.getReadName()); // sanity check
        Assert.assertEquals(163, expectedFor14.getFlags()); // sanity check
        expectedFor14.setFlags(99);
        expectedFor14.setReferenceIndex(5);
        expectedFor14.setAlignmentStart(31356204 + 1);
        expectedFor14.setMappingQuality(19);
        expectedFor14.setCigarString("151M");
        expectedFor14.setMateReferenceIndex(5);
        expectedFor14.setMateAlignmentStart(31356606);
        expectedFor14.setInferredInsertSize(549);
        checkOutputContains(expectedFor14, results);

        SAMRecord expectedFor15 = records.get(15).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", expectedFor15.getReadName()); // sanity check
        Assert.assertEquals(83, expectedFor15.getFlags()); // sanity check
        expectedFor15.setFlags(147);
        expectedFor15.setReferenceIndex(5);
        expectedFor15.setAlignmentStart(31356606);
        expectedFor15.setMappingQuality(0);
        expectedFor15.setCigarString("48M3I100M");
        expectedFor15.setMateReferenceIndex(5);
        expectedFor15.setMateAlignmentStart(31356204 + 1);
        expectedFor15.setInferredInsertSize(-549);
        checkOutputContains(expectedFor15, results);

        // input[16]: A00624:8:HHKYHDSXX:4:2543:5737:19288, reverse strand
        // BWA returns 121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0
        SAMRecord expectedFor16 = records.get(16).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", expectedFor16.getReadName()); // sanity check
        Assert.assertEquals(121, expectedFor16.getFlags()); // sanity check
        expectedFor16.setFlags(121);
        expectedFor16.setReferenceIndex(5);
        expectedFor16.setAlignmentStart(29943927 + 1);
        expectedFor16.setMappingQuality(60);
        expectedFor16.setCigarString("151M");
        expectedFor16.setMateReferenceIndex(5);
        expectedFor16.setMateAlignmentStart(29943927 + 1);
        expectedFor16.setInferredInsertSize(0);
        SAMRecord results16 = results.stream()
                .filter(samRecord -> samRecord.getReadName().equals(expectedFor16.getReadName()))
                .filter(samRecord -> samRecord.getFlags() == 121)
                .findFirst().orElseThrow();
        check(expectedFor16, results16);

        // input[17]: A00624:8:HHKYHDSXX:4:2543:5737:19288, forward strand
        // BWA returns 181,-1,-1,-1,-1,-1,0,0,0,0,"",null,null,5,29943927,0 - not mapped
        SAMRecord expectedFor17 = records.get(17).deepCopy();
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", expectedFor17.getReadName()); // sanity check
        Assert.assertEquals(181, expectedFor17.getFlags()); // sanity check
        expectedFor17.setFlags(181);
        expectedFor17.setReferenceIndex(5);
        expectedFor17.setAlignmentStart(29943927 + 1);
        expectedFor17.setMappingQuality(0);
        expectedFor17.setCigarString("*");
        expectedFor17.setMateReferenceIndex(5);
        expectedFor17.setMateAlignmentStart(29943927 + 1);
        expectedFor17.setInferredInsertSize(0);
        // We have to search for this record manually as it has the same start location as that of its mate.
        SAMRecord results17 = results.stream()
                .filter(samRecord -> samRecord.getReadName().equals(expectedFor17.getReadName()))
                .filter(samRecord -> samRecord.getFlags() == 181)
                .findFirst().orElseThrow();
        check(expectedFor17, results17);
    }

    private void checkOutputContains(SAMRecord expected, List<SAMRecord> output)
    {
        check(expected, findMatchingRecord(expected, output));
    }

    private SAMRecord findMatchingRecord(SAMRecord toFind, List<SAMRecord> output)
    {
        Optional<SAMRecord> result = output.stream()
                .filter(record ->
                        toFind.getReadName().equals(record.getReadName()) && toFind.getAlignmentStart() == record.getAlignmentStart())
                .findFirst();
        if(result.isEmpty())
        {
            Assert.fail("Could not find matching record for " + toFind);
        }
        return result.get();
    }
}

class DummyConfig extends AltContigRemapperConfig
{
    private final PairAligner aligner;

    DummyConfig(PairAligner aligner, String inputBam, String outputBam)
    {
        super(inputBam, outputBam, null);
        this.aligner = aligner;
    }

    @Override
    public PairAligner aligner()
    {
        return aligner;
    }
}

/**
 * This class implements Aligner by hard-coding the alignments for the hla alt records
 * in the test data file tiny.sam to the values returned by a BwaAligner instance that uses
 * the hg38_no_alts reference sequence.
 */
class HlaAligner implements PairAligner
{

    private final List<SAMRecord> records;
    private final Map<Integer, List<BwaMemAlignment>> recordIndexToAlignments = new HashMap<>();

    HlaAligner(List<SAMRecord> records)
    {
        this.records = records;

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(4).getReadName()); // Sanity check
        // 81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0
        // 161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0
        recordIndexToAlignments.put(4, List.of(bwa("81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0")));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(7).getReadName()); // Sanity check
        recordIndexToAlignments.put(7, List.of(bwa("161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0")));

        // 81,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,5,31354347,-564
        // 161,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,5,31354760,564
        // 2225,1,32916241,32916273,42,74,13,0,32,28,42S32M77S,32,null,5,31354760,0
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(8).getReadName()); // Sanity check
        Assert.assertTrue(records.get(8).getReadNegativeStrandFlag()); // Sanity check
        recordIndexToAlignments.put(8, List.of(
                bwa("81,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,5,31354347,-564")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(9).getReadName()); // Sanity check
        Assert.assertFalse(records.get(9).getReadNegativeStrandFlag()); // Sanity check
        recordIndexToAlignments.put(9, List.of(
                bwa("161,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,5,31354760,564"),
                bwa("2225,1,32916241,32916273,42,74,13,0,32,28,42S32M77S,32,null,5,31354760,0")
        ));

        // 97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719
        // 145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(10).getReadName()); // Sanity check
        Assert.assertFalse(records.get(10).getReadNegativeStrandFlag()); // Sanity check
        recordIndexToAlignments.put(10, List.of(
                bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(11).getReadName()); // Sanity check
        recordIndexToAlignments.put(11, List.of(
                bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719")
        ));

        // 97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289
        // 2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0
        // 145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(12).getReadName()); // Sanity check
        recordIndexToAlignments.put(12, List.of(
                bwa("145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(13).getReadName()); // Sanity check
        recordIndexToAlignments.put(13, List.of(
                bwa("97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289"),
                bwa("2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", records.get(14).getReadName()); // Sanity check
        recordIndexToAlignments.put(14, List.of(
                bwa("97,5,31356204,31356355,0,151,19,4,131,121,151M,21T0C19C79G28", "chr6,+31271111,151M,6;", "5,31271513,-84542")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", records.get(15).getReadName()); // Sanity check
        recordIndexToAlignments.put(15, List.of(
                bwa("145,5,31271513,31271664,0,151,18,6,121,111,151M,5G8C0T28G42T51C11", "chr6,-31356606,48M3I100M,10;", "5,31356204,84542")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", records.get(16).getReadName()); // Sanity check
        recordIndexToAlignments.put(16, List.of(
                bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", records.get(17).getReadName()); // Sanity check
        recordIndexToAlignments.put(17, List.of(
                bwa("181,-1,-1,-1,-1,-1,0,0,0,0,\"\",null,null,5,29943927,0")
        ));
    }

    public List<BwaMemAlignment> alignSequence(byte[] bases)
    {
        for(int index = 0; index < records.size(); index++)
        {
            final byte[] recordBases = records.get(index).getReadBases();
            if(Arrays.equals(bases, recordBases))
            {
                return recordIndexToAlignments.get(index);
            }
            if(Arrays.equals(bases, Nucleotides.reverseComplementBases(recordBases)))
            {
                return recordIndexToAlignments.get(index);
            }
        }
        throw new IllegalArgumentException("No alignment found for " + new String(bases));
    }

    @Override
    public ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2)
    {
        List<BwaMemAlignment> resultLeft = alignSequence(bases1);
        List<BwaMemAlignment> resultRight = alignSequence(bases2);
        return new ImmutablePair<>(resultLeft, resultRight);
    }
}
