package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;

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
            MateCigar = extractMateCigar(read);
        }
        else
        {
            LowerCigar = extractMateCigar(read);
            MateCigar = readCigar;
        }
    }

    private static String extractMateCigar(final SAMRecord read)
    {
        String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        return mateCigar != null ? mateCigar : NO_CIGAR;
    }

    public String toString() { return format("cigar(%s) mateCigar(%s) readId(%s)", LowerCigar, MateCigar, Read.getReadName()); }
}
