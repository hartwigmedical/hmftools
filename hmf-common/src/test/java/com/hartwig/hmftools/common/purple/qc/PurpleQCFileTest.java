package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Random;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCFileTest {

    @Test
    public void testReadWrite() throws IOException {
        final PurpleQC expected = create();
        assertEquals(expected, PurpleQCFile.fromLines(PurpleQCFile.toLines(expected)));
    }

    @NotNull
    private static PurpleQC create() {
        final Random random = new Random();
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(90)
                .ploidy(3.17)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.FEMALE)
                .deletedGenes(120)
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .build();
    }
}
