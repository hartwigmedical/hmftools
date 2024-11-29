package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
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
        // using the HlaAligner test aligner (rather than one that uses the entire genome).
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

//        Assert.assertEquals(records.size(), results.size());

        // The first four records in our test data file should not get altered.
        Assert.assertEquals(records.get(0), results.get(0));
        Assert.assertEquals(records.get(1), results.get(1));
        Assert.assertEquals(records.get(2), results.get(2));
        Assert.assertEquals(records.get(3), results.get(3));
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
        super(inputBam, outputBam, "/Users/timlavers/work/apps/bamtools/bin/bamtools");
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
