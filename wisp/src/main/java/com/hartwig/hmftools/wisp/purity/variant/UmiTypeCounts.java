package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT_ORIG;

public class UmiTypeCounts
{
    public int RefNone;
    public int RefSingle;
    public int RefDual;
    public int AlleleNone;
    public int AlleleSingle;
    public int AlleleDual;
    public int Other; // non-alt, non-ref eg other alts

    public static final UmiTypeCounts NO_UMI_COUNTS = new UmiTypeCounts(
            0, 0, 0, 0, 0, 0, 0);

    public static UmiTypeCounts fromAttribute(final Object umiTypeCountsRaw)
    {
        if(umiTypeCountsRaw == null)
            return NO_UMI_COUNTS;

        String[] umiValues = ((String)umiTypeCountsRaw).split(",", UMI_TYPE_COUNT);

        if(umiValues.length < UMI_TYPE_COUNT_ORIG) // temporary allowance for the initial 6-way split without an 'Other' count
            return NO_UMI_COUNTS;

        int index = 0;
        return new UmiTypeCounts(
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]),
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index]),
                index < umiValues.length ? Integer.parseInt(umiValues[index]) : 0);
    }

    public UmiTypeCounts(
            final int refNone, final int refSingle, final int refDual, final int alleleNone, final int alleleSingle, final int alleleDual,
            final int other)
    {
        RefNone = refNone;
        RefSingle = refSingle;
        RefDual = refDual;
        AlleleNone = alleleNone;
        AlleleSingle = alleleSingle;
        AlleleDual = alleleDual;
        Other = other;
    }

    public void add(final UmiTypeCounts other)
    {
        RefNone += other.RefNone;
        RefSingle += other.RefSingle;
        RefDual += other.RefDual;
        AlleleNone += other.AlleleNone;
        AlleleSingle += other.AlleleSingle;
        AlleleDual += other.AlleleDual;
        Other += other.Other;
    }

    public UmiTypeCounts()
    {
        this(0, 0, 0, 0, 0, 0, 0);
    }

    public int refTotal() { return RefNone + RefSingle + RefDual; }
    public int alleleTotal() { return AlleleNone + AlleleSingle + AlleleDual; }
    public int total() { return refTotal() + alleleTotal() + Other; }
    public int dualTotal() { return RefDual + AlleleDual; }
}
