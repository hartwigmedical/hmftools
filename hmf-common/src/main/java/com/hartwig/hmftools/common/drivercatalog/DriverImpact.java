package com.hartwig.hmftools.common.drivercatalog;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public enum DriverImpact {
    MISSENSE,
    NONSENSE,
    SPLICE,
    INFRAME,
    FRAMESHIFT;

    public static DriverImpact select(@NotNull SomaticVariant variant) {
        if (isFrameshift(variant)) {
            return FRAMESHIFT;
        }
        if (isNonsense(variant)) {
            return NONSENSE;
        }
        if (isMissense(variant)) {
            return MISSENSE;
        }

        if (isSplice(variant)) {
            return SPLICE;
        }

        if (isInframe(variant)) {
            return INFRAME;
        }

        throw new UnsupportedOperationException();
    }

    static boolean isFrameshift(SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    static boolean isNonsense(SomaticVariant variant) {
        return variant.type() != VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.NONSENSE_OR_FRAMESHIFT;
    }

    static boolean isMissense(SomaticVariant variant) {
        return variant.type() != VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE;
    }

    static boolean isInframe(SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE;
    }

    static boolean isSplice(SomaticVariant variant) {
        return variant.canonicalCodingEffect() == CodingEffect.SPLICE;
    }
}
