package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VirusLikelihoodTypeTest {

    @Test
    public void canExtractVirusDriverLikelihood() {
        List<AnnotatedVirus> annotatedVirus = Lists.newArrayList();

        annotatedVirus.add(ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(Strings.EMPTY)
                .interpretation(Strings.EMPTY)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH).build());

        annotatedVirus.add(ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(Strings.EMPTY)
                .interpretation(Strings.EMPTY)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.LOW).build());

        annotatedVirus.add(ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(Strings.EMPTY)
                .interpretation(Strings.EMPTY)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.UNKNOWN).build());

        assertEquals(annotatedVirus.get(0).virusDriverLikelihoodType(), VirusLikelihoodType.HIGH);
        assertEquals(annotatedVirus.get(1).virusDriverLikelihoodType(), VirusLikelihoodType.LOW);
        assertEquals(annotatedVirus.get(2).virusDriverLikelihoodType(), VirusLikelihoodType.UNKNOWN);

    }

}