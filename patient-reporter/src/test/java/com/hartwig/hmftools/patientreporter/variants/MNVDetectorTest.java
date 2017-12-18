package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantImpl;

import org.junit.Test;

public class MNVDetectorTest {

    @Test
    public void canDetectPotentialMNVs() {
        final SomaticVariant potentialMNV = new SomaticVariantImpl.Builder().chromosome("X").position(10).ref("A").alt("C").build();
        final SomaticVariant nonMNV = new SomaticVariantImpl.Builder().chromosome("Y").position(10).ref("G").alt("T").build();

        final List<SomaticVariant> allVariants =
                Lists.newArrayList(new SomaticVariantImpl.Builder().chromosome("1").position(10).ref("T").alt("C").build(),
                        new SomaticVariantImpl.Builder().chromosome("X").position(9).ref("A").alt("G").build(), potentialMNV,
                        new SomaticVariantImpl.Builder().chromosome("4").position(5).ref("A").alt("T").build(),
                        new SomaticVariantImpl.Builder().chromosome("3").position(5).ref("T").alt("A").build(), nonMNV);

        final List<SomaticVariant> reportedVariants = Lists.newArrayList(potentialMNV, nonMNV);

        final List<SomaticVariant> potentialMNVs = MNVDetector.locatePotentialMNVs(allVariants, reportedVariants);

        assertTrue(potentialMNVs.contains(potentialMNV));
        assertFalse(potentialMNVs.contains(nonMNV));
        assertEquals(1, potentialMNVs.size());

        // KODU: Make sure the reported variants are not changed by the detector!
        assertEquals(2, reportedVariants.size());
    }
}