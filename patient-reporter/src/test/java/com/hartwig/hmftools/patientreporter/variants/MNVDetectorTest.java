package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.junit.Test;

public class MNVDetectorTest {

    @Test
    public void canDetectPotentialMNVs() {
        final SomaticVariant potentialMNV = new SomaticVariant.Builder().chromosome("X").position(10).build();
        final SomaticVariant nonMNV = new SomaticVariant.Builder().chromosome("Y").position(10).build();

        final List<SomaticVariant> allVariants = Lists.newArrayList(
                new SomaticVariant.Builder().chromosome("1").position(10).build(),
                new SomaticVariant.Builder().chromosome("X").position(9).build(), potentialMNV,
                new SomaticVariant.Builder().chromosome("4").position(5).build(),
                new SomaticVariant.Builder().chromosome("3").position(5).build(), nonMNV);

        final List<SomaticVariant> reportedVariants = Lists.newArrayList(potentialMNV, nonMNV);

        final List<SomaticVariant> potentialMNVs = MNVDetector.locatePotentialMNVs(allVariants, reportedVariants);

        assertTrue(potentialMNVs.contains(potentialMNV));
        assertFalse(potentialMNVs.contains(nonMNV));

        // KODU: Make sure the reported variants are not changed by the detector!
        assertEquals(2, reportedVariants.size());
    }
}