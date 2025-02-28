package com.hartwig.hmftools.common.bam;

import static java.lang.String.format;

import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BamReadLite
{
    public final String Id;
    public final int AlignmentStart;
    public final String Cigar;
    public final short Flags;
    public final String Contig;
    public final short MapQuality;
    public final String MateContig;
    public final int MateStart;
    public final byte[] Bases;
    public final byte[] Quals;
    public final int InsertSize;

    // further ideas:
    // - read group can also be stored once by the owner rather than per record
    // - these tags are common and could be stored as shorts: NM:i:0  AS:i:89 XS:i:89
    // - cigar could be null whenever read is fully aligned (most common scenario)
    // - quals could be stored as an array of value:count eg 10x25, 2x37 etc

    public final List<SAMRecord.SAMTagAndValue> Attributes;

    public BamReadLite(final SAMRecord record, boolean skipReadId)
    {
        Id = skipReadId ? null : record.getReadName();
        AlignmentStart = record.getAlignmentStart();
        Cigar = record.getCigarString();
        Flags = (short)record.getFlags();
        Contig = record.getReferenceName();
        MapQuality = (short)record.getMappingQuality();
        MateContig = record.getMateReferenceName();
        MateStart = record.getMateAlignmentStart();
        Bases = record.getReadBases();
        Quals = record.getBaseQualities();
        InsertSize = record.getInferredInsertSize();
        Attributes = record.getAttributes();
    }

    public static SAMRecord from(final BamReadLite readLite, final SAMFileHeader samFileHeader, final String readId)
    {
        SAMRecord record = new SAMRecord(samFileHeader);

        record.setReadName(readLite.Id != null ? readLite.Id : readId);
        record.setReadBases(readLite.Bases);
        record.setBaseQualities(readLite.Quals);
        record.setMappingQuality(readLite.MapQuality);
        record.setReferenceName(readLite.Contig);
        record.setAlignmentStart(readLite.AlignmentStart);
        record.setCigarString(readLite.Cigar);

        record.setFlags(readLite.Flags);

        record.setMateReferenceName(readLite.MateContig);
        record.setMateAlignmentStart(readLite.MateStart);
        record.setInferredInsertSize(readLite.InsertSize);

        for(SAMRecord.SAMTagAndValue attribute : readLite.Attributes)
        {
            record.setAttribute(attribute.tag, attribute.value);
        }

        return record;
    }

    public String toString()
    {
        return format("id(%s) coords(%s:%d) cigar(%s) mate(%s:%d) flags(%d) insertSize(%d) attributes(%d)",
                Id, Contig, AlignmentStart, Cigar, MateContig, MateStart, Flags, InsertSize, Attributes.size());
    }
}
