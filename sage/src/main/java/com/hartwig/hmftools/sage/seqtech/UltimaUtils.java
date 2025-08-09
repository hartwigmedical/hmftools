package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.CYCLE_BASES;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;

import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public final class UltimaUtils
{
    public static final HashMap<Byte, Integer> CYCLE_BASE_INDEX;

    // ULTIMA CHECK: Access modifiers (?)
    public static final byte INVALID_BASE = -1;

    static
    {
        CYCLE_BASE_INDEX = Maps.newHashMap();
        for(int i = 0; i < CYCLE_BASES.length; i++)
        {
            CYCLE_BASE_INDEX.put(CYCLE_BASES[i], i);
        }
    }

    public static List<Integer> coreHomopolymerLengths(final VariantReadContext readContext)
    {
        List<Homopolymer> homopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        return homopolymers.stream().map(x -> x.Length).collect(Collectors.toList());
    }

    @VisibleForTesting
    public static boolean isCleanSnv(final VariantReadContext readContext)
    {
        if(!readContext.variant().isSNV())
        {
            return false;
        }

        if(readContext.coreLength() != readContext.refBasesBytes().length)
        {
            return false;
        }

        for(int i = readContext.CoreIndexStart; i <= readContext.CoreIndexEnd; ++i)
        {
            if(i == readContext.VarIndex)
            {
                continue;  // the SNV base
            }

            byte readBase = readContext.readBasesBytes()[i];
            byte refBase = readContext.refBasesBytes()[i - readContext.CoreIndexStart];
            if(readBase != refBase)
            {
                return false;
            }
        }

        return true;
    }
}
