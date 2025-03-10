package com.hartwig.hmftools.common.bam;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.READ_GROUP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;

import java.util.List;

import com.google.common.collect.Lists;

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

    // common attribute values
    public final short AlignmentScore;
    public final short NmCount;
    public final short XScore;

    public final SAMRecord.SAMTagAndValue[] Attributes;

    // optimisations
    // read and read group ID may be stored by the caller to avoid duplication across fragments
    // CIGAR elements and alignment blocks are discarded
    // common attributes are extracted from the attributes list and stored as primitives

    // further ideas:
    // - cigar could be null whenever read is fully aligned (most common scenario)
    // - quals could be stored as an array of value:count eg 10x25, 2x37 etc

    public BamReadLite(final SAMRecord record, boolean skipReadIds)
    {
        Id = skipReadIds ? null : record.getReadName();
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

        // convert common attributes to primitives
        List<SAMRecord.SAMTagAndValue> otherAttributes = null;
        short nmCount = -1;
        short alignmentScore = -1;
        short xScore = -1;

        for(SAMRecord.SAMTagAndValue attribute : record.getAttributes())
        {
            if(attribute.tag.equals(ALIGNMENT_SCORE_ATTRIBUTE))
            {
                alignmentScore = intAttributeToShort(attribute);
            }
            else if(attribute.tag.equals(NUM_MUTATONS_ATTRIBUTE))
            {
                nmCount = intAttributeToShort(attribute);
            }
            else if(attribute.tag.equals(XS_ATTRIBUTE))
            {
                xScore = intAttributeToShort(attribute);
            }
            else if(attribute.tag.equals(READ_GROUP_ATTRIBUTE) && skipReadIds)
            {
                continue;
            }
            else
            {
                if(otherAttributes == null)
                    otherAttributes = Lists.newArrayListWithCapacity(3);

                otherAttributes.add(attribute);
            }
        }

        if(otherAttributes != null)
        {
            Attributes = new SAMRecord.SAMTagAndValue[otherAttributes.size()];

            for(int i = 0; i < otherAttributes.size(); ++i)
            {
                Attributes[i] = otherAttributes.get(i);
            }
        }
        else
        {
            Attributes = null;
        }

        NmCount = nmCount;
        AlignmentScore = alignmentScore;
        XScore = xScore;
    }

    private static short intAttributeToShort(final SAMRecord.SAMTagAndValue attribute)
    {
        return ((Integer)attribute.value).shortValue();
    }

    public static SAMRecord from(final BamReadLite readLite, final SAMFileHeader samFileHeader, final String readId, final String rgId)
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

        if(readLite.Attributes != null)
        {
            for(int i = 0; i < readLite.Attributes.length; ++i)
            {
                record.setAttribute(readLite.Attributes[i].tag, readLite.Attributes[i].value);
            }
        }

        if(rgId != null && !rgId.isEmpty())
            record.setAttribute(READ_GROUP_ATTRIBUTE, rgId);

        if(readLite.NmCount >= 0)
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, (int)readLite.NmCount);

        if(readLite.AlignmentScore >= 0)
            record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, (int)readLite.AlignmentScore);

        if(readLite.XScore >= 0)
            record.setAttribute(XS_ATTRIBUTE, (int)readLite.XScore);

        return record;
    }

    public String toString()
    {
        return format("id(%s) coords(%s:%d) cigar(%s) mate(%s:%d) flags(%d) insertSize(%d) attributes(%d)",
                Id, Contig, AlignmentStart, Cigar, MateContig, MateStart, Flags, InsertSize, Attributes != null ? Attributes.length : 0);
    }
}
