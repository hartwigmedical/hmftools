package com.hartwig.hmftools.chord.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

public class SvContext
{
    public final StructuralVariantType mType;
    @Nullable public final Integer mLength;
    public final Integer mLowerInterval;
    public final Integer mUpperInterval;

    public static final List<StructuralVariantType> SV_TYPES = List.of(DEL, DUP, INV, BND);
    private static final List<StructuralVariantType> SV_TYPES_WITH_LENGTH = List.of(DEL, DUP, INV);

    private static final DecimalFormat SV_LENGTH_FORMAT = new DecimalFormat("0E00", new DecimalFormatSymbols(Locale.ENGLISH));
    private static final int[] SV_LENGTH_INTERVALS = { 0, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, Integer.MAX_VALUE};

    private static final String TRANSLOCATION = "TRA";

    private SvContext(
            StructuralVariantType type,
            @Nullable Integer length,
            int lowerInterval,
            int upperInterval
    )
    {
        mType = type;
        mLength = length;
        mLowerInterval = lowerInterval;
        mUpperInterval = upperInterval;
    }

    public static SvContext from(StructuralVariant sv)
    {
        int lowerInterval = Integer.MIN_VALUE;
        int upperInterval = Integer.MIN_VALUE;

        if(sv.Type != BND)
        {
            for(int i = 0; i < SV_LENGTH_INTERVALS.length; i++)
            {
                upperInterval = SV_LENGTH_INTERVALS[i];
                if(sv.Length < upperInterval)
                {
                    lowerInterval = SV_LENGTH_INTERVALS[i-1];
                    break;
                }
            }
        }

        return new SvContext(sv.Type, sv.Length, lowerInterval, upperInterval);
    }

    public String getContextName()
    {
        if(mType == BND)
            return TRANSLOCATION;

        String lowerIntervalString = SV_LENGTH_FORMAT.format(mLowerInterval);
        String upperIntervalString = (mUpperInterval != Integer.MAX_VALUE) ? SV_LENGTH_FORMAT.format(mUpperInterval) : "Inf";

        // Use lower case 'e' for exponent
        lowerIntervalString = lowerIntervalString.replace('E', 'e');
        upperIntervalString = upperIntervalString.replace('E', 'e');

        return String.join("_", mType.toString(), lowerIntervalString, upperIntervalString, "bp");
    }

    public static Map<String, Integer> initializeCounts()
    {
        Map<String, Integer> bins = new LinkedHashMap<>();

        for(StructuralVariantType type : SV_TYPES_WITH_LENGTH)
        {
            for(int i = 0; i < SV_LENGTH_INTERVALS.length-1; i++)
            {
                int lowerInterval = SV_LENGTH_INTERVALS[i];
                int upperInterval = SV_LENGTH_INTERVALS[i + 1];

                String contextName = new SvContext(type, null, lowerInterval, upperInterval).getContextName();

                bins.put(contextName, 0);
            }
        }

        bins.put(TRANSLOCATION, 0);

        return bins;
    }
}
