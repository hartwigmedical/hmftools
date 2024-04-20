package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;
import java.util.Objects;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class SupplementaryReadData
{
    public final String Chromosome;
    public final int Position;
    public final char Strand;
    public final String Cigar;
    public final int MapQuality;
    public final int NM;

    public static final String ALIGNMENTS_DELIM = ";";

    private static final String SUPP_DELIM = ",";
    private static final int SUPP_FIELD_COUNT = 6;

    public static final char SUPP_POS_STRAND = '+';
    public static final char SUPP_NEG_STRAND = '-';

    public SupplementaryReadData(final String chromosome, final int position, final char strand, final String cigar, final int mapQuality,
            final int nm)
    {
        Chromosome = chromosome;
        Position = position;
        Strand = strand;
        Cigar = cigar;
        MapQuality = mapQuality;
        NM = nm;
    }

    public SupplementaryReadData(final String chromosome, final int position, final char strand, final String cigar, final int mapQuality)
    {
        this(chromosome, position, strand, cigar, mapQuality, 0);
    }

    public byte orientation() { return Strand == SUPP_POS_STRAND ? POS_ORIENT : NEG_ORIENT; }

    @Nullable
    @VisibleForTesting
    public static SupplementaryReadData fromAlignment(final String alignment, final String delimiter)
    {
        String[] items = alignment.split(delimiter);
        if(items == null || items.length != SUPP_FIELD_COUNT)
        {
            return null;
        }

        return new SupplementaryReadData(items[0], Integer.parseInt(items[1]), items[2].charAt(0), items[3], Integer.parseInt(items[4]), Integer.parseInt(items[5]));
    }

    @Nullable
    @VisibleForTesting
    public static SupplementaryReadData fromAlignment(final String alignment)
    {
        return fromAlignment(alignment, SUPP_DELIM);
    }

    @Nullable
    public static List<SupplementaryReadData> extractAlignments(final SAMRecord record)
    {
        final String alignmentStr = record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE);
        return alignmentStr != null ? extractAlignments(alignmentStr) : null;
    }

    @Nullable
    public static List<SupplementaryReadData> extractAlignments(@Nullable final String suppData)
    {
        if(suppData == null || suppData.isEmpty())
            return null;

        // example data: 2,33141317,+,94S57M,5,0;
        // but also be multiple: 7,152184341,-,23S32M1I41M54S,0,6;11,66229611,+,115S32M4S,0,0;
        // chr6,6068632,-,35M108S,0,0;chr3,5435688,-,23S39M81S,0,1;chr3,136963678,-,101S31M11S,0,0;

        if(suppData.contains(SUPP_DELIM))
        {
            final List<SupplementaryReadData> output = Lists.newArrayList();

            // do not discard empty strings at the end when splitting
            final String[] alignments = suppData.split(ALIGNMENTS_DELIM, -1);
            final int endIndex = suppData.charAt(suppData.length() - 1) == ALIGNMENTS_DELIM.charAt(0) ?
                    alignments.length - 2 : alignments.length - 1;

            for(int i = 0; i <= endIndex; ++i)
            {
                final String alignment = alignments[i];
                final SupplementaryReadData suppReadData = fromAlignment(alignment);

                if(suppReadData == null)
                    return null;

                output.add(suppReadData);
            }

            return output;
        }
        else
        {
            // alternative delimitation only use for testing
            return List.of(fromAlignment(suppData, ALIGNMENTS_DELIM));
        }
    }

    @Nullable
    public static SupplementaryReadData extractAlignment(final SAMRecord record)
    {
        final String alignmentStr = record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE);
        return alignmentStr != null ? extractAlignment(alignmentStr) : null;
    }

    @Nullable
    public static SupplementaryReadData extractAlignment(@Nullable final String suppData)
    {
        if(suppData == null || suppData.isEmpty())
            return null;

        // example data: 2,33141317,+,94S57M,5,0;
        // but also be multiple: 7,152184341,-,23S32M1I41M54S,0,6;11,66229611,+,115S32M4S,0,0;
        // chr6,6068632,-,35M108S,0,0;chr3,5435688,-,23S39M81S,0,1;chr3,136963678,-,101S31M11S,0,0;

        if(suppData.contains(SUPP_DELIM))
        {
            final String[] alignments = suppData.split(ALIGNMENTS_DELIM, 2);
            return fromAlignment(alignments[0]);
        }
        else
        {
            // Alternative delimitation only use for testing.
            return fromAlignment(suppData, ALIGNMENTS_DELIM);
        }
    }

    public static int alignmentCount(final SAMRecord record)
    {
        if(!record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return 0;

        return alignmentCount(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    public static int alignmentCount(@Nullable final String suppData)
    {
        return suppData != null ? (suppData.contains(SUPP_DELIM) ? suppData.split(ALIGNMENTS_DELIM).length : 1) : 0;
    }

    public static String alignmentsToSamTag(final List<SupplementaryReadData> alignments)
    {
        StringJoiner sj = new StringJoiner(ALIGNMENTS_DELIM);
        alignments.forEach(x -> sj.add(x.asSamTag()));
        return sj.toString();
    }

    public String asCsv()
    {
        return String.format("%s;%d;%c;%s;%d;%d", Chromosome, Position, Strand, Cigar, MapQuality, NM);
    }

    public String asSamTag()
    {
        return String.format("%s,%d,%c,%s,%d,%d", Chromosome, Position, Strand, Cigar, MapQuality, NM);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof SupplementaryReadData))
        {
            return false;
        }
        final SupplementaryReadData that = (SupplementaryReadData) o;
        return Position == that.Position && Strand == that.Strand && MapQuality == that.MapQuality && NM == that.NM
                && Objects.equals(Chromosome, that.Chromosome) && Objects.equals(Cigar, that.Cigar);
    }

    public String toString()
    {
        return String.format("location(%s:%d) strand(%c) cigar(%s) mq(%d) nm(%d)",
                Chromosome, Position, Strand, Cigar, MapQuality, NM);
    }
}
