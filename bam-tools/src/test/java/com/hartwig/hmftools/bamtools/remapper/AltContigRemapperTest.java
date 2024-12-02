package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.*;

import static java.lang.Integer.parseInt;

public class AltContigRemapperTest {

    private final LinkedList<SAMRecord> records = readTestFile();

    @Test
    public void isHlaAltReferenceTest() {
        Assert.assertTrue(AltContigRemapper.isHlaAltReference("HLA-A*01:03"));
        Assert.assertTrue(AltContigRemapper.isHlaAltReference("hla-whatever"));
        Assert.assertFalse(AltContigRemapper.isHlaAltReference("hla whatever"));
        Assert.assertFalse(AltContigRemapper.isHlaAltReference("chrUn_KI270590v1"));
    }

    @Test
    public void hasAltReference() {
        Assert.assertEquals(14, records.size()); // Sanity check

        Assert.assertFalse(AltContigRemapper.hasAltReference(records.get(0)));
        Assert.assertFalse(AltContigRemapper.hasAltReference(records.get(1)));
        Assert.assertTrue(AltContigRemapper.hasAltReference(records.get(7)));
        Assert.assertTrue(AltContigRemapper.hasAltReference(records.get(8)));
        Assert.assertTrue(AltContigRemapper.hasAltReference(records.get(11)));
    }

    @Test
    public void mateHasAltReference() {
        Assert.assertFalse(AltContigRemapper.mateHasAltReference(records.get(0)));
        Assert.assertTrue(AltContigRemapper.mateHasAltReference(records.get(4)));
        Assert.assertTrue(AltContigRemapper.mateHasAltReference(records.get(5)));
        Assert.assertFalse(AltContigRemapper.mateHasAltReference(records.get(7)));
        Assert.assertTrue(AltContigRemapper.mateHasAltReference(records.get(8)));
        Assert.assertTrue(AltContigRemapper.mateHasAltReference(records.get(10)));
    }

    @Test
    public void createRemappedRecordTest() {
        var record = records.get(12);
        // A00624:8:HHKYHDSXX:4:1304:7075:14998	163	HLA-B*08:19N	2778	0	151M	=	3023	289	TGGAAAAGGAGG...C	FFFFF...F	NM:i:0	MD:Z:151	MC:Z:107S44M	MQ:i:0	AS:i:151	XS:i:151

        // This has two alignments in the hg38_no_alts genome:
        // samFlag, refId, refStart, refEnd,   seqStart, seqEnd, mapQual, nMismatches, alignerScore,  suboptimalScore, cigar,     mdTag, xaTag, mateRefId, mateRefStart, templateLen
        // 16,      5,     31354375, 31354419, 0,        44,     60,      1,           39,            20,              44M107S,   22C21,  null, -1,        -1,           0
        // 2048,    1,     32916297, 32916333, 67,       103,    0,       0,           36,            36,              67S36M48S, 36,     null, -1,        -1,           0

        // The first of these gets turned into the following AlignData object:
        // chr6:31354376-31354419	ref location
        // 0-44 raw sequence start and end
        // 0-43 sequence start and end
        // 60	map quality
        // 44M107S	cigar
        // -1	orientation
        // 44	aligned bases
        // 39	score
        // 16	flags
        // 1	nMatches
        // 22C21	xaTag
        // 44	mdTag
        // 0	false
//        var alignData = new AlignData(new ChrBaseRegion("chr6", 31354376, 31354419), 0, 44,
//                60, 39, 16, "44M107S", 1, "22C21", "44");
        BwaMemAlignment alignment = new BwaMemAlignment(16,5,31354375,31354419,0, 44, 60, 1, 39, 20, "44M107S", "22C1", null, -1, -1, 0);
        SAMRecord remappedRecord = AltContigRemapper.createRemappedRecord(record, alignment);


        Assert.assertEquals(record.getReadName(), remappedRecord.getReadName());
        Assert.assertArrayEquals(record.getReadBases(), remappedRecord.getReadBases());
        Assert.assertArrayEquals(record.getBaseQualities(), remappedRecord.getBaseQualities());
        Assert.assertEquals("chr6", remappedRecord.getReferenceName());
        Assert.assertEquals(alignment.getRefStart(), remappedRecord.getAlignmentStart());
//        Assert.assertEquals(alignment.getRefEnd(), remappedRecord.getAlignmentEnd());
        Assert.assertEquals(alignment.getMapQual(), remappedRecord.getMappingQuality());
        Assert.assertEquals(-1, remappedRecord.getMateReferenceIndex().intValue());
        Assert.assertEquals(-1, remappedRecord.getMateAlignmentStart());
    }

    @Test
    public void runTest() {
        // Create a remapper that will remap the tiny bam test data file
        // using the HlaAligner (see below) test aligner (rather than one that uses the entire genome).
        File tempDir = new File("/Users/timlavers/work/junk");
        File outputFile = new File(tempDir, "test.bam");
        File inpuFile = getTestFile("tiny.bam");
        Aligner aligner = new HlaAligner(records);
        AltContigRemapperConfig config = new DummyConfig(aligner, inpuFile.getAbsolutePath(), outputFile.getAbsolutePath());
        AltContigRemapper remapper = new AltContigRemapper(config);

        // Run the alignment.
        remapper.run();

        // Read in the input file as a list of SAM records.
        SamReader samReader = SamReaderFactory.makeDefault().open(outputFile);
        List<SAMRecord> results = new ArrayList<>();
        samReader.iterator().forEachRemaining(results::add);

        // Two of the alignment sequences returned have 2 elements,
        // so there are two more lines in the output than in the input.
        Assert.assertEquals(records.size() + 2, results.size());

        // The first four records in our test data file should not get altered and are already in order.
        check(records.get(0), results.get(0));
        check(records.get(1), results.get(1));
        check(records.get(2), results.get(2));
        check(records.get(3), results.get(3));

        // Expected[4] is input[7] with main alignment remapped.
        // Main alignment: HLA-A*01:03, 3188 -> chr1, 194358862 (see See A00624:8:HHKYHDSXX:2:1516:12156:4225)
        // Mate alignment unchanged.
        SAMRecord expected4 = records.get(7).deepCopy();
        expected4.setReferenceIndex(0);
        expected4.setAlignmentStart(194358862);
        expected4.setMappingQuality(60);
        expected4.setCigarString("35S116M");
        check(expected4, results.get(4));

        // Expected[5] is input[4] with mate remapped
        // Main alignment: unchanged
        // Mate alignment: HLA-A*01:03, 3188 -> chr1, 194358862 (see A00624:8:HHKYHDSXX:2:1516:12156:4225)
        SAMRecord expected5 = records.get(4).deepCopy();
        expected5.setMateReferenceIndex(0);
        expected5.setMateAlignmentStart(194358862);
        check(expected5, results.get(5));

        // Expected[6] is input[8] remapped
        // Main alignment: HLA-B*27:05:18, 3006 -> chr2, 32916241 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*40:06:01:01, 2533 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected6 = records.get(8).deepCopy();
        expected6.setReferenceIndex(1);
        expected6.setAlignmentStart(32916241);
        expected6.setMappingQuality(16);
        expected6.setCigarString("42S32M77S");
        expected6.setMateReferenceIndex(5);
        expected6.setMateAlignmentStart(31354760);
        check(expected6, results.get(6));

        // Expected[7] is input[5] remapped
        // Main alignment: unchanged
        // Mate alignment: HLA-B*40:06:01:01, 2533 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected7 = records.get(5).deepCopy();
        expected7.setMateReferenceIndex(5);
        expected7.setMateAlignmentStart(31354760);
        check(expected7, results.get(7));

        // Expected[8] is input[13] remapped
        // Both the main alignment and the mate refer to HLA-B*08:19N
        // which gets remapped to chr2 for the main alignment and to
        // chr6 for the mate.
        // See the mappings for A00624:8:HHKYHDSXX:4:1304:7075:14998 below.
        SAMRecord expected8 = records.get(13).deepCopy();
        expected8.setReferenceIndex(1);
        expected8.setAlignmentStart(32916297);
        expected8.setMappingQuality(0);
        expected8.setCigarString("67S36M48S");
        expected8.setMateReferenceIndex(5);
        expected8.setMateAlignmentStart(31354513);
        check(expected8, results.get(8));

        // Expected[9] is a second remapping of input[8]
        // Main alignment: HLA-B*27:05:18 -> chr6, 31354347 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*40:06:01:01 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected9 = records.get(8).deepCopy();
        expected9.setReferenceIndex(5);
        expected9.setAlignmentStart(31354347);
        expected9.setMappingQuality(60);
        expected9.setCigarString("72M79S");
        expected9.setMateReferenceIndex(5);
        expected9.setMateAlignmentStart(31354760);
        check(expected9, results.get(9));

        // Expected[10] is a second remapping of input[13]
        // Main alignment: HLA-B*08:19N -> chr6, 31354375 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        // Mate alignment: HLA-B*08:19N -> chr6, 31354513 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        SAMRecord expected10 = records.get(13).deepCopy();
        expected10.setReferenceIndex(5);
        expected10.setAlignmentStart(31354375);
        expected10.setMappingQuality(60);
        expected10.setCigarString("44M107S");
        expected10.setMateReferenceIndex(5);
        expected10.setMateAlignmentStart(31354513);
        check(expected10, results.get(10));

        // Element 12 is the result of remapping record 13.
        // Main alignment: HLA-B*08:19N -> chr6, 31354513 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        // Mate alignment: HLA-B*08:19N -> chr6, 31354513 (see A00624:8:HHKYHDSXX:4:1304:7075:14998)
        SAMRecord expected11 = records.get(12).deepCopy();
        expected11.setReferenceIndex(5);
        expected11.setAlignmentStart(31354513);
        expected11.setMappingQuality(60);
        expected11.setCigarString("151M");
        expected11.setMateReferenceIndex(5);
        expected11.setMateAlignmentStart(31354375);
        check(expected11, results.get(11));

        // Element 13 is the result of remapping record 10.
        // Main alignment: HLA-B*40:06:01:01 -> chr6, 31354760 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        // Mate alignment: HLA-B*27:05:18 -> chr6, 31354347 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected12 = records.get(9).deepCopy();
        expected12.setReferenceIndex(5);
        expected12.setAlignmentStart(31354760);
        expected12.setMappingQuality(60);
        expected12.setCigarString("151M");
        expected12.setMateReferenceIndex(5);
        expected12.setMateAlignmentStart(31354347);
        check(expected12, results.get(12));

        // Element 14 is the result of remapping record 11.
        // Main alignment: HLA-B*07:05:01 -> chr6, 31355729 (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        // Mate alignment: HLA-B*07:05:01 -> chr6, 31356297 (see A00624:8:HHKYHDSXX:1:2559:3224:21292)
        SAMRecord expected13 = records.get(11).deepCopy();
        expected13.setReferenceIndex(5);
        expected13.setAlignmentStart(31355729);
        expected13.setMappingQuality(60);
        expected13.setCigarString("151M");
        expected13.setMateReferenceIndex(5);
        expected13.setMateAlignmentStart(31356297);
        check(expected13, results.get(13));

        // Element 15 is the result of remapping record 10.
        // Main alignment: HLA-B*07:05:01, 711 -> chr6, 31356297  (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        // Mate alignment: HLA-B*07:05:01, 1279 -> chr6, 31355729  (see A00624:8:HHKYHDSXX:1:1446:18213:29684)
        SAMRecord expected14 = records.get(10).deepCopy();
        expected14.setReferenceIndex(5);
        expected14.setAlignmentStart(31356297);
        expected14.setMappingQuality(60);
        expected14.setCigarString("151M");
        expected14.setMateReferenceIndex(5);
        expected14.setMateAlignmentStart(31355729);
        check(expected14, results.get(14));

        // Element 16 is the result of remapping record 6.
        // Main alignment: chr6_GL000250v2_alt is unchanged
        // Mate alignment: HLA-A*01:03, 3188 -> chr1, 194358862  (see A00624:8:HHKYHDSXX:2:1516:12156:4225)
        SAMRecord expected15 = records.get(6).deepCopy();
        expected15.setMateReferenceIndex(0);
        expected15.setMateAlignmentStart(194358862);
        check(expected15, results.get(15));
    }

    private void check(SAMRecord expected, SAMRecord actual) {
        // We do the same checks as would be done by using SAMRecord.equals()
        // except that we do not check attributes.
        Assert.assertEquals(expected.getAlignmentStart(), actual.getAlignmentStart());
        Assert.assertEquals(expected.getFlags(), actual.getFlags());
        Assert.assertEquals(expected.getInferredInsertSize(), actual.getInferredInsertSize());
        Assert.assertEquals(expected.getMappingQuality(), actual.getMappingQuality());
        Assert.assertEquals(expected.getMateAlignmentStart(), actual.getMateAlignmentStart());
        Assert.assertEquals(expected.getMateReferenceIndex(), actual.getMateReferenceIndex());
        Assert.assertEquals(expected.getReferenceIndex(), actual.getReferenceIndex());
        Assert.assertEquals(expected.getReadName(), actual.getReadName());
//        Assert.assertEquals(expected.getAttributes(), actual.getAttributes());
        Assert.assertArrayEquals(expected.getBaseQualities(), actual.getBaseQualities());
        Assert.assertEquals(expected.getCigar(), actual.getCigar());
        Assert.assertEquals(expected.getMateReferenceName(), actual.getMateReferenceName());
        Assert.assertArrayEquals(expected.getReadBases(), actual.getReadBases());
        Assert.assertEquals(expected.getReferenceName(), actual.getReferenceName());
    }

    private LinkedList<SAMRecord> readTestFile() {
        var samReader = SamReaderFactory. makeDefault().open(getTestFile("tiny.sam"));
        var records = new LinkedList<SAMRecord>();
        samReader.forEach(r -> records.add(r));
        return records;
    }

    private File getTestFile(String name) {
        ClassLoader classLoader = getClass().getClassLoader();
        var fileName = classLoader.getResource(name).getFile();
        return new File(fileName);
    }
}

class DummyConfig extends AltContigRemapperConfig {
    private Aligner aligner;

    DummyConfig(Aligner aligner, String inputBam, String outputBam) {
        super(inputBam, outputBam, "/Users/timlavers/work/apps/samtools-1.21/samtools");
        this.aligner = aligner;
    }

    @Override
    public Aligner aligner() {
        return aligner;
    }
}


/**
 * This class implements Aligner by hard-coding the alignments for the hla alt records
 * in the test data file tiny.sam to the values returned by a BwaAligner instance that uses
 * the hg38_no_alts reference sequence.
 */
class HlaAligner implements Aligner {

    private final List<SAMRecord> records;
    private Map<Integer, List<BwaMemAlignment>> recordIndexToAlignments = new HashMap<>();

    HlaAligner(List<SAMRecord> records) {
        this.records = records;

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(7).getReadName()); // Sanity check
        recordIndexToAlignments.put(7, List.of(bwa("0,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,-1,-1,0")));

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
    public List<BwaMemAlignment> alignSequence(byte[] bases) {
        for (int index = 0; index < records.size(); index++) {
            if (Arrays.equals(bases, records.get(index).getReadBases())) {
                return recordIndexToAlignments.get(index);
            }
        }
        throw new IllegalArgumentException("No alignment found for " + Arrays.toString(bases));
    }

    private BwaMemAlignment bwa(String data) {
        String[] parts = data.split(",");
        return new BwaMemAlignment(
                parseInt(parts[0]),
                parseInt(parts[1]),
                parseInt(parts[2]),
                parseInt(parts[3]),
                parseInt(parts[4]),
                parseInt(parts[5]),
                parseInt(parts[6]),
                parseInt(parts[7]),
                parseInt(parts[8]),
                parseInt(parts[9]),
                parts[10],
                parts[11],
                parts[12],
                parseInt(parts[13]),
                parseInt(parts[14]),
                parseInt(parts[15])
                );
    }
}
