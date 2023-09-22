package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;

public class UmiTypeCounts
{
    public int RefNone;
    public int RefSingle;
    public int RefDual;
    public int AlleleNone;
    public int AlleleSingle;
    public int AlleleDual;

    public static final UmiTypeCounts NO_UMI_COUNTS = new UmiTypeCounts(
            0, 0, 0, 0, 0, 0);

    public static UmiTypeCounts fromAttribute(final Object umiTypeCountsRaw)
    {
        if(umiTypeCountsRaw == null)
            return NO_UMI_COUNTS;

        String[] umiValues = ((String)umiTypeCountsRaw).split(",", UMI_TYPE_COUNT);

        if(umiValues.length != UMI_TYPE_COUNT)
            return NO_UMI_COUNTS;

        int index = 0;
        return new UmiTypeCounts(
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]),
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index]));
    }

    public UmiTypeCounts(
            final int refNone, final int refSingle, final int refDual, final int alleleNone, final int alleleSingle, final int alleleDual)
    {
        RefNone = refNone;
        RefSingle = refSingle;
        RefDual = refDual;
        AlleleNone = alleleNone;
        AlleleSingle = alleleSingle;
        AlleleDual = alleleDual;
    }

    public void add(final UmiTypeCounts other)
    {
        RefNone += other.RefNone;
        RefSingle += other.RefSingle;
        RefDual += other.RefDual;
        AlleleNone += other.AlleleNone;
        AlleleSingle += other.AlleleSingle;
        AlleleDual += other.AlleleDual;
    }

    public UmiTypeCounts()
    {
        this(0, 0, 0, 0, 0, 0);
    }

    public int refTotal() { return RefNone + RefSingle + RefDual; }
    public int alleleTotal() { return AlleleNone + AlleleSingle + AlleleDual; }
    public int total() { return refTotal() + alleleTotal(); }
    public int dualTotal() { return RefDual + AlleleDual; }
}
