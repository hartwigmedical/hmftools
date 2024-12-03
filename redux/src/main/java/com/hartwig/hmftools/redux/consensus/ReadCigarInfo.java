package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadCigarInfo
{
    public final SAMRecord Read;
    public final String LowerCigar;
    public final String MateCigar;

    public ReadCigarInfo(final SAMRecord read, boolean isLowerRead, boolean isPrimary)
    {
        Read = read;

        String readCigar;

        if(isPrimary)
        {
            readCigar = read.getCigarString();
        }
        else
        {
            readCigar = SupplementaryReadData.extractAlignment(read).Cigar;
        }

        if(isLowerRead)
        {
            LowerCigar = readCigar;
            MateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        }
        else
        {
            LowerCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            MateCigar = readCigar;
        }
    }

    public String toString() { return format("cigar(%s) mateCigar(%s) readId(%s)", LowerCigar, MateCigar, Read.getReadName()); }
}
