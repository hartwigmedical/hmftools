package com.hartwig.hmftools.common.bam.testutilities;

import java.util.Arrays;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class PairedRecordsBuilder
{
    private final byte BASE_QUAL_F = (byte) ('A' + 5);
    private final String readName;
    private final SAMFileHeader header;

    public PairedRecordsBuilder(final String readName, final SAMFileHeader header)
    {
        this.readName = readName;
        this.header = header;
    }

    public Pair<SAMRecord, SAMRecord> build(Pair<BasesRegion, BasesRegion> baseRegionPair)
    {
        BasesRegion readRegion = baseRegionPair.getLeft();
        BasesRegion mateRegion = baseRegionPair.getRight();
        int length = readRegion.mBases.length;
        SAMRecord read1 = createRecord(length, true);
        // flags: 1 (paired), 2, (read mapped in proper pair), 32 (mate reverse strand), 64 (1st in pair)
        read1.setFlags(99);
        int chrIndex = header.getSequenceIndex(baseRegionPair.getLeft().chromosome());
        int mateChrIndex = header.getSequenceIndex(baseRegionPair.getRight().chromosome());
        read1.setReferenceIndex(chrIndex);
        read1.setAlignmentStart(readRegion.start());
        read1.setReadBases(readRegion.mBases);
        read1.setMateAlignmentStart(mateRegion.start());
        read1.setMateReferenceIndex(mateChrIndex);

        SAMRecord read2 = createRecord(length, false);
        // flags 1, 2, 16 (read reverse strand), 128 (2nd in pair)
        read2.setFlags(147);
        read2.setReferenceIndex(mateChrIndex);
        read2.setAlignmentStart(mateRegion.start());
        read2.setReadBases(mateRegion.mBases);
        read2.setMateAlignmentStart(readRegion.start());
        read2.setMateReferenceIndex(chrIndex);

        return Pair.of(read1, read2);
    }

    private SAMRecord createRecord(int length, boolean forward)
    {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(readName);
        record.setHeader(header);
        record.setCigarString(length + "M");
        record.setMappingQuality(60);
        record.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, 0);
        record.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, "0");
        record.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, length);
        record.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, length + "M");
        int lengthSign = forward ? 1 : -1;
        record.setInferredInsertSize(lengthSign * 2 * length);
        record.setAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE, 60); // 60?
        // Perfect base qualities.
        record.setBaseQualities(baseQualities(length));

        return record;
    }

    private byte[] baseQualities(int length)
    {
        byte[] qualities = new byte[length];
        Arrays.fill(qualities, BASE_QUAL_F);
        return qualities;
    }
}
