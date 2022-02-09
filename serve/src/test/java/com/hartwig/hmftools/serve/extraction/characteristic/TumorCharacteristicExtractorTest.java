package com.hartwig.hmftools.serve.extraction.characteristic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.util.Strings;
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
    public void canDetermineCutOffMSI() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.tumorCharacteristicAnnotation());
        assertEquals("MSI >= 4", characteristic.cutoff());

    }

    @Test
    public void canDetermineCutOffMSS() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.tumorCharacteristicAnnotation());
        assertEquals("MSS < 4", characteristic.cutoff());
    }

    @Test
    public void canDetermineCutOffTMLLow() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
        assertEquals("TML < 140", characteristic.cutoff());
    }

    @Test
    public void canDetermineCutOffTMLHigh() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
        assertEquals("TML >= 140", characteristic.cutoff());
    }

    @Test
    public void canDetermineCutOffHRD() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.tumorCharacteristicAnnotation());
        assertEquals("HRD >= 0.5", characteristic.cutoff());
    }

    @Test
    public void canDetermineCutOffSource() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB, "TML >= 200");
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
        assertEquals("TML >= 200", characteristic.cutoff());
    }

    @Test
    public void canExtractMicrosatelliteUnstableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractMicrosatelliteStableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractHighTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractLowTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TMB, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractHrDeficientCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractHPVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HPV, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HPV_POSITIVE, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractEBVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, EBV, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.EBV_POSITIVE, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canExtractImmunoHlaCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, IMMUNO_HLA, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.IMMUNO_HLA, characteristic.tumorCharacteristicAnnotation());
    }

    @Test
    public void canFilterUnknownCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, "Not a tumor characteristic", Strings.EMPTY));
    }

    @Test
    public void canFilterWrongTypes() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.COMPLEX, MSI, Strings.EMPTY));
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