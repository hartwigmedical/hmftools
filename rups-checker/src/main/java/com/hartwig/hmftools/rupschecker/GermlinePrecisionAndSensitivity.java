package com.hartwig.hmftools.rupschecker;

import com.hartwig.hmftools.common.variant.VariantType;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import org.jetbrains.annotations.NotNull;

public class GermlinePrecisionAndSensitivity extends PrecisionAndSensitivity {
    
    @NotNull
    private final String config;
    
    public GermlinePrecisionAndSensitivity(VariantType variantType, String algorithm, String config) {
        super(variantType, algorithm);
        this.config = config;
    }
    
    public String getConfig() {
        return config;
    }    
    
    @Override
    public String getHeader(boolean produceMDOutput) {
        String delimit = delimiter(produceMDOutput);
        StringBuilder builder = new StringBuilder();
        builder.append(delimit).append("Config").
                append(super.getHeader(produceMDOutput));   
        return builder.toString();
    }
    
    @Override
    public String getValueString(boolean produceMDOutput) {
        String delimit = delimiter(produceMDOutput);
        StringBuilder builder = new StringBuilder();
        builder.append(delimit).append(config).
                append(super.getHeader(produceMDOutput));   
        return builder.toString();
    }
}
