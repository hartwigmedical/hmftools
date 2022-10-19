package com.hartwig.hmftools.serve.extraction.characteristic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristic;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristicsComparator;

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
    public void canDetermineMSI() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, MSI);

        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.name());
        assertNull(characteristic.comparator());
        assertNull(characteristic.cutoff());
    }

    @Test
    public void canDetermineMSS() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, MSS);

        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.name());
        assertNull(characteristic.comparator());
        assertNull(characteristic.cutoff());
    }

    @Test
    public void canDetermineTMLLow() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, LOW_TML + " < 140");

        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.LOWER, characteristic.comparator());
        assertEquals(140, characteristic.cutoff(), EPSILON);
    }

    @Test
    public void canDetermineTMLHigh() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, HIGH_TML + " >= 140");

        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_GREATER, characteristic.comparator());
        assertEquals(140, characteristic.cutoff(), EPSILON);
    }

    @Test
    public void canDetermineTMBLow() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, LOW_TMB + " <= 3");

        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.EQUAL_OR_LOWER, characteristic.comparator());
        assertEquals(3, characteristic.cutoff(), EPSILON);
    }

    @Test
    public void canDetermineTMBHigh() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, HIGH_TMB + " > 14.5");

        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, characteristic.name());
        assertEquals(TumorCharacteristicsComparator.GREATER, characteristic.comparator());
        assertEquals(14.5, characteristic.cutoff(), EPSILON);
    }

    @Test
    public void canDetermineHRD() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, HRD);

        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.name());
        assertNull(characteristic.comparator());
        assertNull(characteristic.cutoff());
    }

    @Test
    public void canExtractHPVPositiveCharacteristic() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, HPV);

        assertEquals(TumorCharacteristicAnnotation.HPV_POSITIVE, characteristic.name());
        assertNull(characteristic.comparator());
        assertNull(characteristic.cutoff());
    }

    @Test
    public void canExtractEBVPositiveCharacteristic() {
        TumorCharacteristicExtractor extractor = buildTestExtractor();

        TumorCharacteristic characteristic = extractor.extract(EventType.CHARACTERISTIC, EBV);

        assertEquals(TumorCharacteristicAnnotation.EBV_POSITIVE, characteristic.name());
        assertNull(characteristic.comparator());
        assertNull(characteristic.cutoff());
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
                Sets.newHashSet(HIGH_TML),
                Sets.newHashSet(LOW_TML),
                Sets.newHashSet(HIGH_TMB),
                Sets.newHashSet(LOW_TMB),
                Sets.newHashSet(HRD),
                Sets.newHashSet(HPV),
                Sets.newHashSet(EBV));
    }
}