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
    private static final String MSS = "mss";
    private static final String HIGH_TMB = "high_tmb";
    private static final String LOW_TMB = "low_tmb";
    private static final String HRD = "hrd";
    private static final String HPV = "hpv";
    private static final String EBV = "ebv";
    private static final String IMMUNO_HLA = "hla";

    @Test
    public void canExtractMicrosatelliteUnstableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.MICROSATELLITE_UNSTABLE, characteristic);
    }

    @Test
    public void canExtractMicrosatelliteStableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.MICROSATELLITE_STABLE, characteristic);
    }

    @Test
    public void canExtractHighTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic);
    }

    @Test
    public void canExtractLowTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TMB);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.LOW_TUMOR_MUTATIONAL_LOAD, characteristic);
    }

    @Test
    public void canExtractHrDeficientCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic);
    }

    @Test
    public void canExtractHPVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HPV);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.HPV_POSITIVE, characteristic);
    }

    @Test
    public void canExtractEBVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, EBV);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.EBV_POSITIVE, characteristic);
    }

    @Test
    public void canExtractImmunoHlaCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, IMMUNO_HLA);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristic.IMMUNO_HLA, characteristic);
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
        return new TumorCharacteristicExtractor(Sets.newHashSet(MSI),
                Sets.newHashSet(MSS),
                Sets.newHashSet(HIGH_TMB),
                Sets.newHashSet(LOW_TMB),
                Sets.newHashSet(HRD),
                Sets.newHashSet(HPV),
                Sets.newHashSet(EBV),
                Sets.newHashSet(IMMUNO_HLA));
    }
}