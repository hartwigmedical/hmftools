package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.common.variant.SageVcfTags.CONSENSUS_TAG_TYPE_COUNT;

public class UmiTypeCounts
{
    public int TotalNone;
    public int TotalSingle;
    public int TotalDual;
    public int AlleleNone;
    public int AlleleSingle;
    public int AlleleDual;

    public static final UmiTypeCounts NO_UMI_COUNTS = new UmiTypeCounts(
            0, 0, 0, 0, 0, 0);

    public static UmiTypeCounts fromAttribute(final Object umiTypeCountsRaw)
    {
        if(umiTypeCountsRaw == null)
            return NO_UMI_COUNTS;

        String[] umiValues = ((String)umiTypeCountsRaw).split(",", CONSENSUS_TAG_TYPE_COUNT);

        if(umiValues.length < CONSENSUS_TAG_TYPE_COUNT)
            return NO_UMI_COUNTS;

        int index = 0;
        return new UmiTypeCounts(
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]),
                Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index++]), Integer.parseInt(umiValues[index]));
    }

    public UmiTypeCounts(
            final int totalNone, final int totalSingle, final int totalDual, final int alleleNone, final int alleleSingle, final int alleleDual)
    {
        TotalNone = totalNone;
        TotalSingle = totalSingle;
        TotalDual = totalDual;
        AlleleNone = alleleNone;
        AlleleSingle = alleleSingle;
        AlleleDual = alleleDual;
    }

    public void add(final UmiTypeCounts other)
    {
        TotalNone += other.TotalNone;
        TotalSingle += other.TotalSingle;
        TotalDual += other.TotalDual;
        AlleleNone += other.AlleleNone;
        AlleleSingle += other.AlleleSingle;
        AlleleDual += other.AlleleDual;
    }

    public UmiTypeCounts()
    {
        this(0, 0, 0, 0, 0, 0);
    }

    public int alleleCount() { return AlleleNone + AlleleSingle + AlleleDual; }
    public int totalCount() { return TotalNone + TotalSingle + TotalDual; }
}
