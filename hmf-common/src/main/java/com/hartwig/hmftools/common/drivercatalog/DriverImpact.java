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

    @NotNull
    public static DriverImpact select(@NotNull SomaticVariant variant) {
        if (isFrameshift(variant)) {
            return FRAMESHIFT;
        } else if (isNonsense(variant)) {
            return NONSENSE;
        } else if (isMissense(variant)) {
            return MISSENSE;
        } else if (isSplice(variant)) {
            return SPLICE;
        } else if (isInframe(variant)) {
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
