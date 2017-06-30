package com.hartwig.hmftools.rupschecker;

import com.hartwig.hmftools.common.variant.VariantType;

public class SomaticPrecisionAndSensitivity extends PrecisionAndSensitivity {
    
    public SomaticPrecisionAndSensitivity(VariantType variantType, String algorithm) {
        super(variantType, algorithm);
    }
}
