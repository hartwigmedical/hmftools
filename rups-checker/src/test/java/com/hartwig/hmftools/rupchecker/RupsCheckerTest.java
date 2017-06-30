package com.hartwig.hmftools.rupchecker;


import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.rupschecker.PrecisionAndSensitivity;
import com.hartwig.hmftools.rupschecker.SomaticPrecisionAndSensitivity;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class RupsCheckerTest {

    private static final double EPSILON = 1.0e-4;
    
    @Test
    public void calculateSensitivity() {
        
        // sensitivity = TP / (TP + FN)
        
        PrecisionAndSensitivity stats = new SomaticPrecisionAndSensitivity(VariantType.INDEL, "foobar");
        
        stats.incTruePositive();
        assertEquals(1, stats.TP());
        
        stats.incFalseNegatives();
        assertEquals(1, stats.FN());
        assertEquals(0.5, stats.sensitivity(), EPSILON);
        
        stats.incFalseNegatives();
        stats.incFalseNegatives();
        assertEquals(3, stats.FN());
        assertEquals(0.25, stats.sensitivity(), EPSILON);
        
        stats.incFalsePositives();  // NOP to sensi calc
        assertEquals(1, stats.FP());
        assertEquals(0.25, stats.sensitivity(), EPSILON);
    }
    
    @Test
    public void calculatePrecision() {
        
        // precision = TP / (TP + FP)
        
        PrecisionAndSensitivity stats = new SomaticPrecisionAndSensitivity(VariantType.SNP, "foobar");
        
        stats.incTruePositive();
        assertEquals(1, stats.TP());
         
        stats.incFalsePositives();
        assertEquals(1, stats.FP());
        assertEquals(0.5, stats.precision(), EPSILON);
        
        stats.incFalsePositives();
        stats.incFalsePositives();
        assertEquals(3, stats.FP());
        assertEquals(0.25, stats.precision(), EPSILON);
        
        stats.incFalseNegatives();  // NOP to precision calc
        assertEquals(1, stats.FN());
        assertEquals(0.25, stats.precision(), EPSILON);
    }
}
