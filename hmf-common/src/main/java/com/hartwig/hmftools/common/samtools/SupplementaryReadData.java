package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import htsjdk.samtools.SAMRecord;

public class SupplementaryReadData
{
    public final String Chromosome;
    public final int Position;
    public final char Strand;
    public final String Cigar;
    public final int MapQuality;

    private static final String SUPP_DELIM = ",";
    private static final String ALIGNMENTS_DELIM = ";";
    private static final int SUPP_FIELD_COUNT = 6;

    public static final char SUPP_POS_STRAND = '+';
    public static final char SUPP_NEG_STRAND = '-';

    public static SupplementaryReadData from(final SAMRecord record)
    {
        String alignmentStr = record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE);
        return alignmentStr != null ? from(alignmentStr) : null;
    }

    public static SupplementaryReadData from(final String suppData)
    {
        if(suppData == null)
            return null;

        // example data: 2,33141317,+,94S57M,5,0;
        // but also be multiple: 7,152184341,-,23S32M1I41M54S,0,6;11,66229611,+,115S32M4S,0,0;
        // chr6,6068632,-,35M108S,0,0;chr3,5435688,-,23S39M81S,0,1;chr3,136963678,-,101S31M11S,0,0;

        String[] items = null;
        if(suppData.contains(SUPP_DELIM))
        {
            // return the first alignment
            final String[] alignments = suppData.split(ALIGNMENTS_DELIM);
            items = alignments[0].split(SUPP_DELIM);
        }
        else
        {
            // alternative delimination - only used for testing
            items = suppData.split(ALIGNMENTS_DELIM);
        }

        if(items == null || items.length != SUPP_FIELD_COUNT)
            return null;

        return new SupplementaryReadData(items[0], Integer.parseInt(items[1]), items[2].charAt(0), items[3], Integer.parseInt(items[4]));
    }

    public static int alignmentCount(final SAMRecord record)
    {
        if(!record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return 0;

        return alignmentCount(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    public static int alignmentCount(final String suppData)
    {
        return suppData != null ? (suppData.contains(SUPP_DELIM) ? suppData.split(ALIGNMENTS_DELIM).length : 1) : 0;
    }

    public SupplementaryReadData(final String chromosome, final int position, final char strand, final String cigar, final int mapQuality)
    {
        Chromosome = chromosome;
        Position = position;
        Strand = strand;
        Cigar = cigar;
        MapQuality = mapQuality;
    }

    public String toString()
    {
        return String.format("location(%s:%d) strand(%c) cigar(%s) mq(%d)",
            Chromosome, Position, Strand, Cigar, MapQuality);
    }

    public String asCsv() { return String.format("%s;%d;%c;%s;%d;0", Chromosome, Position, Strand, Cigar, MapQuality); }
    public String asSamTag() { return String.format("%s,%d,%c,%s,%d,0", Chromosome, Position, Strand, Cigar, MapQuality); }
}
