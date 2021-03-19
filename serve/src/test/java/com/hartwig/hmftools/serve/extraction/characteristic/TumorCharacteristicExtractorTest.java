package com.hartwig.hmftools.serve.extraction.characteristic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorCharacteristicExtractorTest {

    private static final String MSI = "msi";
    private static final String TMB = "tmb";
    private static final String HRD = "hrd";

    @Test
    public void canExtractMicrosatelliteUnstableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.MICROSATELLITE_UNSTABLE, characteristic);
    }

    @Test
    public void canExtractHighTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, TMB);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic);
    }

    @Test
    public void canExtractHrDeficientCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic);
    }

    @Test
    public void canFilterUnknownCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, "Not a tumor characteristic"));
    }

    @Test
    public void canFilterWrongTypes() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.COMPLEX, MSI));
    }

    @NotNull
    private static TumorCharacteristicExtractor buildTestExtractor() {
        return new TumorCharacteristicExtractor(Sets.newHashSet(MSI), Sets.newHashSet(TMB), Sets.newHashSet(HRD));
    }
}