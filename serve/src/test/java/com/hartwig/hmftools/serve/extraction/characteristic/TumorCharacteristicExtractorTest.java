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

    private static final double EPSILON = 1e-10;

    private static final String MSI = "msi";
    private static final String MSS = "mss";
    private static final String HIGH_TML = "high_tml";
    private static final String LOW_TML = "low_tml";
    private static final String HIGH_TMB = "high_tmb";
    private static final String LOW_TMB = "low_tmb";
    private static final String HRD = "hrd";
    private static final String HPV = "hpv";
    private static final String EBV = "ebv";

    @Test
    public void canDetermineCutOffMSI() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, "MSI >= 4"));
        assertEquals(4,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, "MSI >= 4"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffMSS() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.LOWER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, "MSS < 4"));
        assertEquals(4,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, "MSS < 4"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMLLow() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TML, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.LOWER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, "TML < 140"));
        assertEquals(140,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, "TML < 140"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMLHigh() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 140"));
        assertEquals(140,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 140"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMBLow() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.LOWER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, "TMB < 3"));
        assertEquals(3,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, "TMB < 3"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMBHigh() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, "TMB >= 14"));
        assertEquals(14,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, "TMB >= 14"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffHRD() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        "HRD >= 0.5"));
        assertEquals(0.5,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        "HRD >= 0.5"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffSource() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, "TML >= 200");
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.name());

        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 200"));
        assertEquals(200.0,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 200"),
                EPSILON);

    }

    @Test
    public void canExtractMicrosatelliteUnstableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, null));
        assertEquals(4.0,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, null),
                EPSILON);
    }

    @Test
    public void canExtractMicrosatelliteStableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.LOWER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, null));
        assertEquals(4.0, TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, null), EPSILON);
    }

    @Test
    public void canExtractHighTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, null));
        assertEquals(140,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, null),
                EPSILON);
    }

    @Test
    public void canExtractLowTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TML, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.LOWER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, null));
        assertEquals(140,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, null),
                EPSILON);
    }

    @Test
    public void canExtractHrDeficientCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER,
                TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, null));
        assertEquals(0.5,
                TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, null),
                EPSILON);
    }

    @Test
    public void canExtractHPVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HPV, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HPV_POSITIVE, characteristic.name());
        assertNull(TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HPV_POSITIVE, null));
        assertNull(TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HPV_POSITIVE, null));
    }

    @Test
    public void canExtractEBVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, EBV, null);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.EBV_POSITIVE, characteristic.name());
        assertNull(TumorCharacteristicExtractor.determineComparator(TumorCharacteristicAnnotation.HPV_POSITIVE, null));
        assertNull(TumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HPV_POSITIVE, null));
    }

    @Test
    public void canFilterUnknownCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, "Not a tumor characteristic", null));
    }

    @Test
    public void canFilterWrongTypes() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.COMPLEX, MSI, null));
    }

    @NotNull
    private static TumorCharacteristicExtractor buildTestExtractor() {
        return new TumorCharacteristicExtractor(Sets.newHashSet(MSI),
                Sets.newHashSet(MSS),
                Sets.newHashSet(HIGH_TML),
                Sets.newHashSet(LOW_TML),
                Sets.newHashSet(HIGH_TMB),
                Sets.newHashSet(LOW_TMB),
                Sets.newHashSet(HRD),
                Sets.newHashSet(HPV),
                Sets.newHashSet(EBV));
    }
}