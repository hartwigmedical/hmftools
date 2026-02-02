package com.hartwig.hmftools.amber.contamination;

import static org.junit.Assert.assertEquals;

public class TumorContaminationFileTest
{
    /*
    @Test
    public void testReadWrite()
    {
        TumorContamination contamination = create(new Random());
        List<String> write = TumorContaminationFile.toLines(Lists.newArrayList(contamination));
        List<TumorContamination> read = TumorContaminationFile.fromLines(write);
        assertEquals(contamination, read.get(0));
    }

    @NotNull
    private static TumorContamination create(@NotNull Random random)
    {
        final BaseDepth template = createRandom(random);

        final BaseDepthData normalDepth = ModifiableBaseDepth.create()
                .from(template)
                .setReadDepth(random.nextInt())
                .setRefSupport(random.nextInt())
                .setAltSupport(random.nextInt());

        final BaseDepthData tumorDepth = ModifiableBaseDepth.create()
                .from(template)
                .setReadDepth(random.nextInt())
                .setRefSupport(random.nextInt())
                .setAltSupport(random.nextInt());

        return ImmutableTumorContamination.builder().from(template).normal(normalDepth).tumor(tumorDepth).build();
    }

    static BaseDepth createRandom(@NotNull final Random random)
    {
        return createRandom(String.valueOf(random.nextInt(22)), random);
    }

    static BaseDepth createRandom(@NotNull final String chromosome, @NotNull final Random random)
    {
        return ModifiableBaseDepth.create()
                .setChromosome(chromosome)
                .setPosition(random.nextInt())
                .setRef(BaseDepth.Base.A)
                .setAlt(BaseDepth.Base.T)
                .setReadDepth(random.nextInt())
                .setRefSupport(random.nextInt())
                .setAltSupport(random.nextInt())
                .setIndelCount(0);
    }
     */
}
