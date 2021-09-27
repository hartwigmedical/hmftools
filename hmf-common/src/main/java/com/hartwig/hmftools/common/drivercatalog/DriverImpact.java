package com.hartwig.hmftools.common.drivercatalog;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public enum DriverImpact
{
    MISSENSE,
    NONSENSE,
    SPLICE,
    INFRAME,
    FRAMESHIFT,
    UNKNOWN;

    @NotNull
    public static DriverImpact select(@NotNull SomaticVariant variant)
    {
        return select(variant.type(), variant.canonicalCodingEffect());
    }

    @NotNull
    public static DriverImpact select(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        if(isFrameshift(variantType, canonicalCodingEffect))
        {
            return FRAMESHIFT;
        }
        else if(isNonsense(variantType, canonicalCodingEffect))
        {
            return NONSENSE;
        }
        else if(isMissense(variantType, canonicalCodingEffect))
        {
            return MISSENSE;
        }
        else if(isSplice(canonicalCodingEffect))
        {
            return SPLICE;
        }
        else if(isInframe(variantType, canonicalCodingEffect))
        {
            return INFRAME;
        }

        return UNKNOWN;
    }

    public static boolean isFrameshift(SomaticVariant variant)
    {
        return isFrameshift(variant.type(), variant.canonicalCodingEffect());
    }

    public static boolean isNonsense(SomaticVariant variant)
    {
        return isNonsense(variant.type(), variant.canonicalCodingEffect());
    }

    public static boolean isMissense(SomaticVariant variant)
    {
        return isMissense(variant.type(), variant.canonicalCodingEffect());
    }

    public static boolean isInframe(SomaticVariant variant)
    {
        return isInframe(variant.type(), variant.canonicalCodingEffect());
    }

    public static boolean isSplice(SomaticVariant variant)
    {
        return isSplice(variant.canonicalCodingEffect());
    }

    static boolean isFrameshift(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType == VariantType.INDEL && canonicalCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    static boolean isNonsense(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType != VariantType.INDEL && canonicalCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    static boolean isMissense(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType != VariantType.INDEL && canonicalCodingEffect == CodingEffect.MISSENSE;
    }

    static boolean isInframe(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType == VariantType.INDEL && canonicalCodingEffect == CodingEffect.MISSENSE;
    }

    static boolean isSplice(CodingEffect canonicalCodingEffect)
    {
        return canonicalCodingEffect == CodingEffect.SPLICE;
    }
}
