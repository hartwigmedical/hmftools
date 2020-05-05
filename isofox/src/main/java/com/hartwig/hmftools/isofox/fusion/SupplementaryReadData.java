package com.hartwig.hmftools.isofox.fusion;

public class SupplementaryReadData
{
    public final String Chromosome;
    public final int Position;
    public final char Strand;
    public final String Cigar;
    public final int MapQuality;

    private static final String SUPP_DELIM = ",";
    private static final int SUPP_FIELD_COUNT = 6;

    public static SupplementaryReadData from(final String suppData)
    {
        // example data: 21;42870046;-;46S30M;255;0;
        final String[] items = suppData.split(SUPP_DELIM);
        if(items.length != SUPP_FIELD_COUNT)
            return null;

        return new SupplementaryReadData(items[0], Integer.parseInt(items[1]), items[2].charAt(0), items[3], Integer.parseInt(items[4]));
    }

    public SupplementaryReadData(final String chromosome, final int position, final char strand, final String cigar, final int mapQuality)
    {
        Chromosome = chromosome;
        Position = position;
        Strand = strand;
        Cigar = cigar;
        MapQuality = mapQuality;
    }
}
