package com.hartwig.hmftools.orange.conversion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;

import org.junit.Test;

public class OrangeConversionTest {

    @Test
    public void shouldConvertAnnotatedVirus() {
        assertTrue(OrangeConversion.convert(TestCommonDatamodelFactory.emptyVirusInterpreterData()).allViruses().isEmpty());
        assertEqualsValue(OrangeConversion.convert(TestCommonDatamodelFactory.minimalVirusInterpreterData()),
                TestCommonDatamodelFactory.minimalVirusInterpreterData());
        assertEqualsValue(OrangeConversion.convert(TestCommonDatamodelFactory.exhaustiveVirusInterpreterData()),
                TestCommonDatamodelFactory.exhaustiveVirusInterpreterData());
    }

    private static void assertEqualsValue(VirusInterpreterData converted, com.hartwig.hmftools.common.virus.VirusInterpreterData input) {
        assertEqualsValue(converted.allViruses(), input.allViruses());
        assertEqualsValue(converted.reportableViruses(), input.reportableViruses());
    }

    private static void assertEqualsValue(List<AnnotatedVirus> converted, List<com.hartwig.hmftools.common.virus.AnnotatedVirus> input) {
        assertEquals(converted.size(), input.size());
        for (int i = 0; i < input.size(); i++) {
            assertEqualsValue(converted.get(i), input.get(i));
        }
    }

    private static void assertEqualsValue(AnnotatedVirus converted, com.hartwig.hmftools.common.virus.AnnotatedVirus input) {
        assertEquals(converted.name(), input.name());
        assertEquals(converted.qcStatus().name(), input.qcStatus().name());
        assertEquals(converted.integrations(), input.integrations());
        assertEquals(converted.percentageCovered(), input.percentageCovered(), 0.01);
        assertEquals(converted.meanCoverage(), input.meanCoverage(), 0.01);
        assertEquals(converted.expectedClonalCoverage(), input.expectedClonalCoverage());
        assertEquals(converted.reported(), input.reported());
        assertEquals(converted.virusDriverLikelihoodType().name(), input.virusDriverLikelihoodType().name());
    }
}