package com.hartwig.hmftools.common.drivercatalog;

import java.util.Comparator;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

public class TsgImpactComparator implements Comparator<SomaticVariant> {

    @Override
    public int compare(final SomaticVariant o1, final SomaticVariant o2) {
        int firstWins = -1;
        int secondWins = 1;

        if (o1.type() == o2.type() && o1.canonicalCodingEffect() == o2.canonicalCodingEffect()) {
            return 0;
        }

        if (isFrameshift(o1)) {
            return firstWins;
        } else if (isFrameshift(o2)) {
            return secondWins;
        }

        if (isNonsense(o1)) {
            return firstWins;
        } else if (isNonsense(o2)) {
            return secondWins;
        }

        if (isSplice(o1)) {
            return firstWins;
        } else if (isSplice(o2)) {
            return secondWins;
        }

        if (isMissense(o1)) {
            return firstWins;
        } else if (isMissense(o2)) {
            return secondWins;
        }

        throw new UnsupportedOperationException();
    }

    static boolean isSplice(SomaticVariant variant) {
        return variant.canonicalCodingEffect() == CodingEffect.SPLICE;
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

}