package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.alignmentsToSamTag;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.CHR_4;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.READ_UNMAPPER;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.checkTransformRead;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMRecord;

public class PartitionReaderTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final TestBamWriter mWriter;

    private final PartitionReader mPartitionReader;

    public PartitionReaderTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES_REPEAT_40);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES_REPEAT_40.length());

        ReduxConfig config = new ReduxConfig(
                mRefGenome, false, false, true, new ReadUnmapper(Collections.emptyMap()));

        mWriter = new TestBamWriter(config);

        mPartitionReader = createPartitionRead(config, mWriter);
    }

    @Test
    public void testMultipleSupplementaries()
    {
        // one primary has multiple supplementaries and the other has a different supp
        int readPos = 100;
        int matePos = 200;
        int suppPos = 1800;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        List<SupplementaryReadData> read1Supps = Lists.newArrayList(
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_2, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(SUPPLEMENTARY_ATTRIBUTE, alignmentsToSamTag(read1Supps));
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos + 2, TEST_READ_BASES, "2S98M", CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_2, suppPos + 1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);
        mPartitionReader.flushReadPositions();

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos + 2, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate2);

        mPartitionReader.processRead(mate1);
        mPartitionReader.processRead(mate2);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.processRead(read1Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_2, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        SAMRecord read2Supp1 = createSamRecord(
                read2.getReadName(), CHR_2, suppPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_2, 1001, 2000));

        mPartitionReader.processRead(read1Supp2);
        mPartitionReader.processRead(read2Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp3 = createSamRecord(
                read1.getReadName(), CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_3, 1001, 2000));

        mPartitionReader.processRead(read1Supp3);
        mPartitionReader.postProcessRegion();

        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount()); // since the supplementaries aren't mapped to the same location
    }

    @Test
    public void testInconsistentSupplementaries()
    {
        // supplementaries arriving out of order and mapped to different locations
        int supp2Pos1 = 950;
        int readPos = 1950;
        int supp2Pos2 = 2950;
        int matePos = 3500;
        int supp1Pos1 = 4500;
        int supp1Pos2 = 5500;

        String readId1 = mReadIdGen.nextId();
        String readId2 = mReadIdGen.nextId();

        SAMRecord read2Supp1 = createSamRecord(
                readId2, CHR_1, supp2Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 999));
        mPartitionReader.processRead(read2Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1 = createSamRecord(
                readId1, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp1Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                readId2, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp2Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        int regionStart = 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);

        mPartitionReader.postProcessRegion();

        SAMRecord read2Supp2 = createSamRecord(
                readId2, CHR_1, supp2Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);
        setSecondInPair(read2Supp2);

        regionStart += 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read2Supp2);
        mPartitionReader.postProcessRegion();

        SAMRecord mate1 = createSamRecord(
                readId1, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp1Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                readId2, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp2Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);
        setSecondInPair(mate2);

        regionStart += 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(mate1);
        mPartitionReader.processRead(mate2);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);
        setSecondInPair(read1Supp1);

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        // note last 2 regions are out of order
        regionStart = 5000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1Supp1);
        mPartitionReader.postProcessRegion();

        regionStart = 4000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1Supp2);
        mPartitionReader.postProcessRegion();

        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());
    }

    @Test
    public void testUnmappedMateRead()
    {
        ReduxConfig config = new ReduxConfig(mRefGenome, false, false, false, READ_UNMAPPER);

        TestBamWriter writer = new TestBamWriter(config);

        PartitionReader partitionReader = createPartitionRead(config, writer);

        // mate is unmapped at original location by RegionUnmapper - simulated by changing its coords and setting the unmapped attribute
        // it is then processed twice at the new location
        SAMRecord read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, 500, false,
                false, null,true, TEST_READ_CIGAR);

        // in an unmapping region
        SAMRecord mate = createSamRecord(
                read.getReadName(), CHR_3, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100, true,
                false, null,false, TEST_READ_CIGAR);

        SAMRecord mateOriginal = cloneSamRecord(mate, mate.getReadName());

        // mate processed by RegionUnmapper at original location
        checkTransformRead(mate, mate.getReferenceName());
        assertTrue(mate.getReadUnmappedFlag());
        assertFalse(mate.getMateUnmappedFlag());

        // mate processed again at original location by partition reader

        partitionReader.setupRegion(new ChrBaseRegion(CHR_3, 1, 1000));
        partitionReader.processRead(mateOriginal);

        assertTrue(mateOriginal.getReadUnmappedFlag());

        partitionReader.postProcessRegion();

        // confirm was unmapped again but not written to BAM
        assertEquals(0, mWriter.nonConsensusWriteCount());

        // read and now-unmapped mate processed at new location
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));
        partitionReader.processRead(read);
        partitionReader.processRead(mate);

        assertTrue(read.getMateUnmappedFlag());
        assertTrue(mate.getReadUnmappedFlag());

        partitionReader.postProcessRegion();

        assertEquals(2, writer.nonConsensusWriteCount());
    }

    @Test
    public void testFullyUnmappedReads()
    {
        // scenario 1 - both reads go from mapped to unmapped and are dropped by the RegionUnmapped and the PartitionReader
        ReduxConfig config = new ReduxConfig(mRefGenome, false, false, false, READ_UNMAPPER);

        TestBamWriter writer = new TestBamWriter(config);

        PartitionReader partitionReader = createPartitionRead(config, writer);

        // mate is unmapped at original location by RegionUnmapper - simulated by changing its coords and setting the unmapped attribute
        // it is then processed twice at the new location
        SAMRecord read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_3, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_4, 1000, false,
                false, null,true, TEST_READ_CIGAR);

        // in an unmapping region
        SAMRecord mate = createSamRecord(
                read.getReadName(), CHR_4, 1000, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, 500, true,
                false, null,false, TEST_READ_CIGAR);

        SAMRecord readOriginal = cloneSamRecord(read, read.getReadName());
        SAMRecord mateOriginal = cloneSamRecord(mate, mate.getReadName());

        // reads are processed by RegionUnmapper at original location
        checkTransformRead(read, read.getReferenceName());
        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());

        checkTransformRead(mate, mate.getReferenceName());
        assertTrue(mate.getReadUnmappedFlag());
        assertTrue(mate.getMateUnmappedFlag());

        // reads are processed again at original location by partition reader

        partitionReader.setupRegion(new ChrBaseRegion(CHR_3, 1, 1000));
        partitionReader.processRead(readOriginal);

        assertTrue(readOriginal.getReadUnmappedFlag());
        assertTrue(readOriginal.getMateUnmappedFlag());

        partitionReader.postProcessRegion();

        partitionReader.setupRegion(new ChrBaseRegion(CHR_4, 1, 4000));
        partitionReader.processRead(mateOriginal);

        assertTrue(mateOriginal.getReadUnmappedFlag());
        assertTrue(mateOriginal.getMateUnmappedFlag());

        partitionReader.postProcessRegion();

        // confirm was unmapped again but not written to BAM
        assertEquals(0, mWriter.nonConsensusWriteCount());

        // scenario 2 - one read is already unmapped and its mate also becomes unmapped
        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_3, 500, TEST_READ_BASES, TEST_READ_CIGAR, CHR_4, 1000, false,
                false, null,false, TEST_READ_CIGAR);
        read.setMateUnmappedFlag(true);

        assertTrue(read.getMateUnmappedFlag());

        // in an unmapping region
        mate = createSamRecord(
                read.getReadName(), CHR_3, 500, TEST_READ_BASES, NO_CIGAR, CHR_3, 500, false,
                false, null,false, TEST_READ_CIGAR);
        mate.setReadUnmappedFlag(true);

        readOriginal = cloneSamRecord(read, read.getReadName());
        mateOriginal = cloneSamRecord(mate, mate.getReadName());

        // reads are processed by RegionUnmapper at original location
        checkTransformRead(read, read.getReferenceName());
        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());

        checkTransformRead(mate, mate.getReferenceName());
        assertTrue(mate.getReadUnmappedFlag());
        assertTrue(mate.getMateUnmappedFlag());

        // reads are processed again at original location by partition reader

        partitionReader.setupRegion(new ChrBaseRegion(CHR_3, 1, 1000));
        partitionReader.processRead(readOriginal);
        partitionReader.processRead(mateOriginal);

        assertTrue(readOriginal.getReadUnmappedFlag());
        assertTrue(readOriginal.getMateUnmappedFlag());

        assertTrue(mateOriginal.getReadUnmappedFlag());
        assertTrue(mateOriginal.getMateUnmappedFlag());

        partitionReader.postProcessRegion();

        // confirm was unmapped again but not written to BAM
        assertEquals(0, mWriter.nonConsensusWriteCount());
    }

    @Test
    public void testUltima()
    {
        ReduxConfig config = new ReduxConfig(mRefGenome, false, false, true, new ReadUnmapper(Collections.emptyMap()), ULTIMA);
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 2 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 2 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 2 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 2 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 3 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150 + 3 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));

        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150, true));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 150, true));

        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 175 + 4 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));

        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 200 + 4 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 200 + 4 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 200 + 5 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 200 + 6 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));
        partitionReader.processRead(createUnpairedRecord(CHR_1, 100, 200 + 6 * SINGLE_END_JITTER_COLLAPSE_DISTANCE, false));

        partitionReader.postProcessRegion();

        assertEquals(19, writer.nonConsensusWriteCount());
        assertEquals(6, writer.consensusWriteCount());
    }

    @Test
    public void testIlluminaPCRClusterCount()
    {
        ReduxConfig config = new ReduxConfig(mRefGenome, false, false, true, new ReadUnmapper(Collections.emptyMap()), ILLUMINA);

        Function<String, SAMRecord> getSamRecord =
                readName -> createSamRecord(readName, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 1000, false, false, null, false, TEST_READ_CIGAR);

        // invalid read names
        TestBamWriter writer = new TestBamWriter(config);
        PartitionReader partitionReader = createPartitionRead(config, writer);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        String readName1 = "READ_001";
        String readName2 = "READ_002";

        partitionReader.processRead(getSamRecord.apply(readName1));
        partitionReader.processRead(getSamRecord.apply(readName2));

        partitionReader.postProcessRegion();

        List<SAMRecord> consensusReads = writer.WrittenRecords.stream()
                .filter(read -> read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                .collect(Collectors.toList());

        assertEquals(1, consensusReads.size());

        SAMRecord consensusRead = consensusReads.get(0);
        String crAttributeValue = consensusRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        assertEquals("2;2", crAttributeValue);

        // two optical duplicates
        writer = new TestBamWriter(config);
        partitionReader = createPartitionRead(config, writer);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        readName1 = "A00624:8:HHKYHDSXX:4:2140:1:1";
        readName2 = "A00624:8:HHKYHDSXX:4:2141:1:1";

        partitionReader.processRead(getSamRecord.apply(readName1));
        partitionReader.processRead(getSamRecord.apply(readName2));

        partitionReader.postProcessRegion();

        consensusReads = writer.WrittenRecords.stream().
                filter(read -> read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                .collect(Collectors.toList());

        assertEquals(1, consensusReads.size());

        consensusRead = consensusReads.get(0);
        crAttributeValue = consensusRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        assertEquals("2;2;1", crAttributeValue);

        // two pcr duplicates
        writer = new TestBamWriter(config);
        partitionReader = createPartitionRead(config, writer);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        readName1 = "A00624:8:HHKYHDSXX:4:2140:1:1";
        readName2 = "A00624:8:HHKYHDSXX:5:2140:1:1";

        partitionReader.processRead(getSamRecord.apply(readName1));
        partitionReader.processRead(getSamRecord.apply(readName2));

        partitionReader.postProcessRegion();

        consensusReads = writer.WrittenRecords.stream()
                .filter(read -> read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                .collect(Collectors.toList());

        assertEquals(1, consensusReads.size());

        consensusRead = consensusReads.get(0);
        crAttributeValue = consensusRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        assertEquals("2;2;2", crAttributeValue);

        // three duplicates on different tiles
        writer = new TestBamWriter(config);
        partitionReader = createPartitionRead(config, writer);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        readName1 = "A00624:8:HHKYHDSXX:4:2140:1:1";
        readName2 = "A00624:8:HHKYHDSXX:4:2141:1:1";
        String readName3 = "A00624:8:HHKYHDSXX:4:2142:1:1";

        partitionReader.processRead(getSamRecord.apply(readName1));
        partitionReader.processRead(getSamRecord.apply(readName2));
        partitionReader.processRead(getSamRecord.apply(readName3));

        partitionReader.postProcessRegion();

        consensusReads = writer.WrittenRecords.stream()
                .filter(read -> read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                .collect(Collectors.toList());

        assertEquals(1, consensusReads.size());

        consensusRead = consensusReads.get(0);
        crAttributeValue = consensusRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        assertEquals("3;3;3", crAttributeValue);

        // seven duplicates with five pcr clusters
        writer = new TestBamWriter(config);
        partitionReader = createPartitionRead(config, writer);
        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        readName1 = "A00624:8:HHKYHDSXX:4:2140:1:1";
        readName2 = "A00624:8:HHKYHDSXX:4:2140:1:1";
        readName3 = "A00624:8:HHKYHDSXX:4:2141:1:1";
        String readName4 = "A00624:8:HHKYHDSXX:4:2141:3001:1";
        String readName5 = "A00624:8:HHKYHDSXX:4:2142:1:1";
        String readName6 = "A00624:8:HHKYHDSXX:4:2142:1:1";
        String readName7 = "A00624:8:HHKYHDSXX:4:2142:3001:1";

        partitionReader.processRead(getSamRecord.apply(readName1));
        partitionReader.processRead(getSamRecord.apply(readName2));
        partitionReader.processRead(getSamRecord.apply(readName3));
        partitionReader.processRead(getSamRecord.apply(readName4));
        partitionReader.processRead(getSamRecord.apply(readName5));
        partitionReader.processRead(getSamRecord.apply(readName6));
        partitionReader.processRead(getSamRecord.apply(readName7));

        partitionReader.postProcessRegion();

        consensusReads = writer.WrittenRecords.stream()
                .filter(read -> read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                .collect(Collectors.toList());

        assertEquals(1, consensusReads.size());

        consensusRead = consensusReads.get(0);
        crAttributeValue = consensusRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        assertEquals("7;7;5", crAttributeValue);
    }

    private static SAMRecord createUnpairedRecord(final String chromosome, final int readStart, int readEnd, boolean isReversed)
    {
        return TestUtils.createUnpairedRecord(READ_ID_GEN.nextId(), chromosome, readStart, readEnd, isReversed);
    }
}
