package com.hartwig.hmftools.common.driver;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

public enum DriverImpact
{
    MISSENSE,
    NONSENSE,
    SPLICE,
    INFRAME,
    FRAMESHIFT,
    UNKNOWN;

    public static DriverImpact select(SomaticVariant variant)
    {
        return select(variant.type(), variant.canonicalCodingEffect());
    }

    public static DriverImpact select(final VariantType variantType, final CodingEffect canonicalCodingEffect)
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

    public static boolean isFrameshift(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType == VariantType.INDEL && canonicalCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    public static boolean isNonsense(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType != VariantType.INDEL && canonicalCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    public static boolean isMissense(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType != VariantType.INDEL && canonicalCodingEffect == CodingEffect.MISSENSE;
    }

    public static boolean isInframe(VariantType variantType, CodingEffect canonicalCodingEffect)
    {
        return variantType == VariantType.INDEL && canonicalCodingEffect == CodingEffect.MISSENSE;
    }

    public static boolean isSplice(CodingEffect canonicalCodingEffect)
    {
        return canonicalCodingEffect == CodingEffect.SPLICE;
    }
}
