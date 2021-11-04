package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Random;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQC;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCFileTest {

    @Test
    public void testReadWrite() {
        final Random random = new Random();
        final PurpleQC expected = create(random);
        assertEquals(expected, PurpleQCFile.fromLines(PurpleQCFile.toLines(expected)));
    }

    @NotNull
    public static PurpleQC create(@NotNull final Random random) {
        return ImmutablePurpleQC.builder()
                .copyNumberSegments(random.nextInt())
                .unsupportedCopyNumberSegments(random.nextInt())
                .purity(random.nextInt(5))
                .amberGender(Gender.values()[random.nextInt(Gender.values().length)])
                .cobaltGender(Gender.values()[random.nextInt(Gender.values().length)])
                .deletedGenes(random.nextInt())
                .contamination(random.nextInt(5))
                .method(FittedPurityMethod.values()[random.nextInt(FittedPurityMethod.values().length)])
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .amberMeanDepth(random.nextInt(100))
                .build();
    }
}
