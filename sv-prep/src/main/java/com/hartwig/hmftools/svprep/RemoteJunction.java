package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

public class RemoteJunction
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public int Support;

    public RemoteJunction(final String chromosome, final int position, final byte orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        Support = 0;
    }

    public static RemoteJunction fromSupplementaryData(final SupplementaryReadData suppData)
    {
        byte orientation = suppData.Strand == '+' ? POS_ORIENT : NEG_ORIENT;
        return new RemoteJunction(suppData.Chromosome, suppData.Position, orientation);
    }

    public boolean matches(final RemoteJunction other)
    {
        return Chromosome.equals(other.Chromosome) && Position == other.Position && Orientation == other.Orientation;
    }

    public boolean matches(final String chromosome, final int position, final byte orientation)
    {
        return Chromosome.equals(chromosome) && Position == position && Orientation == orientation;
    }

    public String toString() { return format("loc(%s:%d:%d) reads(%d)", Chromosome, Position, Orientation, Support); }
}
