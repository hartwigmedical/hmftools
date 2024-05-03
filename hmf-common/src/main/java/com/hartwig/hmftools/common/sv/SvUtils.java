package com.hartwig.hmftools.common.sv;

import static java.lang.Math.abs;

import htsjdk.samtools.SAMRecord;

public final class SvUtils
{
    public static final int DEFAULT_DISCORDANT_FRAGMENT_LENGTH = 1000;

    // must match the small deldup threshold in Esvee
    public static final int SMALL_DELDUP_SIZE = 1000;

    public static boolean isDiscordant(final SAMRecord record) { return isDiscordant(record, DEFAULT_DISCORDANT_FRAGMENT_LENGTH); }

    public static boolean isDiscordant(final SAMRecord record, final int discordantPairFragmentLength)
    {
        if(!record.getReadPairedFlag())
            return false;

        if(!record.getReferenceName().equals(record.getMateReferenceName()))
            return true;

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return true;

        int fragmentSize = abs(record.getInferredInsertSize());

        return fragmentSize == 0 || fragmentSize >= discordantPairFragmentLength;
    }

}
