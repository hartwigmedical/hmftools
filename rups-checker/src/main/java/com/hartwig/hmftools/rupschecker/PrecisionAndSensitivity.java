package com.hartwig.hmftools.rupschecker;

import com.hartwig.hmftools.common.variant.VariantType;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Objects;
import org.jetbrains.annotations.NotNull;

public abstract class PrecisionAndSensitivity {
    
    @NotNull
    protected final VariantType variantType;
    
    @NotNull
    protected final String algorithm;
    
    protected long truePositives;
    
    protected long falsePositives;
    
    protected long falseNegatives;
   
    protected static final String ASCII_DEL = "\t\t";
    protected static final String MD_DEL = "|";
    
    public PrecisionAndSensitivity(VariantType variantType, String algorithm) {
        this.variantType = variantType;
        this.algorithm = algorithm;
        this.truePositives = 0;
        this.falsePositives = 0;
        this.falseNegatives = 0;
    }
    
    public long TP() {
        return truePositives;
    }
    
    public long FP() {
        return falsePositives;
    }
    
    public long FN() {
        return falseNegatives;
    }
    
    public double precision() {
        if (truePositives + falsePositives > 0) {
            return (double)truePositives / (double)(truePositives + falsePositives);
        } else {
            //throw new RuntimeException("precision calculation: TP + FP > 0");
            return 0.0;
        }
    }

    public double sensitivity() {
        if (truePositives + falseNegatives > 0) {
            return (double)truePositives / (double)(truePositives + falseNegatives);
        } else {
            //throw new RuntimeException("sensitivity calculation: TP + FN > 0");
            return 0.0;
        }
    }
    
    public void incTruePositive() {
        truePositives++;
    }
    
    public void incFalsePositives() {
        falsePositives++;
    }
    
    public void incFalseNegatives() {
        falseNegatives++;
    }
    
    public VariantType getVariantType() {
        return variantType;
    }
    
    public String getAlgorithm() {
        return algorithm;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 71 * hash + Objects.hashCode(this.variantType);
        hash = 71 * hash + Objects.hashCode(this.algorithm);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final PrecisionAndSensitivity other = (PrecisionAndSensitivity) obj;
        if (!Objects.equals(this.algorithm, other.algorithm)) {
            return false;
        }
        if (this.variantType != other.variantType) {
            return false;
        }
        return true;
    }
    
    public String getHeader(boolean produceMDOutput) {
        String delimit = delimiter(produceMDOutput);
        StringBuilder builder = new StringBuilder();
        builder.append(delimit).append("Type").
                append(delimit).append("Algo").
                append(delimit).append("TP").
                append(delimit).append("FP").
                append(delimit).append("FN").
                append(delimit).append("Prec").
                append(delimit).append("Sens").
                append(delimit);
        return builder.toString();
    }
    
    public String getValueString(boolean produceMDOutput) {
        String delimit = delimiter(produceMDOutput);
        DecimalFormat df = new DecimalFormat("##.#");
        df.setRoundingMode(RoundingMode.UP);
        StringBuilder builder = new StringBuilder();
        builder.append(delimit).append(variantType).
                append(delimit).append(algorithm).
                append(delimit).append(truePositives).
                append(delimit).append(falsePositives).
                append(delimit).append(falseNegatives).
                append(delimit).append(df.format(precision())).
                append(delimit).append(df.format(sensitivity())).
                append(delimit);
        return builder.toString();
    }
    
    protected String delimiter(boolean produceMDOutput) {
        return produceMDOutput ? MD_DEL : ASCII_DEL;
    }
}
