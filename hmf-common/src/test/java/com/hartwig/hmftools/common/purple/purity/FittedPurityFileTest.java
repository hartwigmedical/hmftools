package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.fromLine;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Random;
import java.util.StringJoiner;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityFileTest {

    @Test
    public void testStripDecimalInTumorMutationalLoad() {
        PurityContext victim = createRandomContext(new Random());
        String line = FittedPurityFile.toString(victim);
        String[] array = line.split("\t");
        array[19] = array[19] + ".9892ADV";
        StringJoiner joiner = new StringJoiner("\t");
        for (String s : array) {
            joiner.add(s);
        }

        assertEquals(victim, FittedPurityFile.fromLine(joiner.toString()));
    }

    @Test
    public void testInputAndOutput() {
        final Random random = new Random();
        final PurityContext input = createRandomContext(random);

        final List<String> lines = toLines(input);
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith("purity"));

        final PurityContext output = fromLine(lines.get(1));
        assertEquals(input, output);
    }

    @Test
    public void testCompatibilityWith2_14() throws IOException {
        FittedPurityFile.fromLine(Resources.readLines(Resources.getResource("purple/v2-47.purple.purity"), Charset.defaultCharset())
                .get(1));
    }

    private static PurityContext createRandomContext(@NotNull Random random) {
        return ImmutablePurityContext.builder()
                .version(random.nextInt() + "a")
                .bestFit(createRandomPurity(random))
                .score(createRandomScore(random))
                .gender(Gender.values()[random.nextInt(Gender.values().length)])
                .status(FittedPurityStatus.values()[random.nextInt(FittedPurityStatus.values().length)])
                .polyClonalProportion(nextDouble(random))
                .wholeGenomeDuplication(random.nextBoolean())
                .microsatelliteIndelsPerMb(random.nextDouble())
                .microsatelliteStatus(MicrosatelliteStatus.MSI)
                .tumorMutationalBurdenPerMb(random.nextDouble())
                .tumorMutationalBurdenStatus(TumorMutationalStatus.HIGH)
                .tumorMutationalLoad(random.nextInt(100_000_000))
                .tumorMutationalLoadStatus(TumorMutationalStatus.LOW)
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .addGermlineAberrations(GermlineAberration.values()[random.nextInt(GermlineAberration.values().length)])
                .build();
    }

    @NotNull
    private static FittedPurityScore createRandomScore(@NotNull Random random) {
        return ImmutableFittedPurityScore.builder()
                .minPurity(nextDouble(random))
                .maxPurity(nextDouble(random))
                .minPloidy(nextDouble(random))
                .maxPloidy(nextDouble(random))
                .minDiploidProportion(nextDouble(random))
                .maxDiploidProportion(nextDouble(random))
                .build();
    }

    @NotNull
    private static FittedPurity createRandomPurity(@NotNull Random random) {
        return createRandomPurityBuilder(random).build();
    }

    @NotNull
    static ImmutableFittedPurity.Builder createRandomPurityBuilder(@NotNull Random random) {
        return ImmutableFittedPurity.builder()
                .purity(nextDouble(random))
                .normFactor(nextDouble(random))
                .score(nextDouble(random))
                .diploidProportion(nextDouble(random))
                .ploidy(nextDouble(random))
                .somaticPenalty(nextDouble(random));
    }

    private static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000D) / 10000D;
    }
}
