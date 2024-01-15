package com.hartwig.hmftools.esvee.old;

import java.util.Arrays;

public class PrimaryPhasingTest
{
    private static SupportedAssembly assembly(final String... supportFragments)
    {
        final SupportedAssembly assembly = new SupportedAssembly("Dummy", "ATCG") {
            @Override
            public String toString()
            {
                return Arrays.toString(supportFragments);
            }
        };
        // FIXME: This
        //assembly.SupportFragments.addAll(Arrays.asList(supportFragments));
        return assembly;
    }

    /* CHASHA FIXME
    @Test
    public void canPhaseSimple()
    {
        final SupportedAssembly one = assembly("A", "B", "C");
        final SupportedAssembly two = assembly("D", "E", "F");
        final SupportedAssembly three = assembly("A", "B", "C");
        final SupportedAssembly four = assembly("H");

        final List<Set<SupportedAssembly>> sets = PrimaryPhasing.run(List.of(one, two, three, four));
        assertTrue(sets).hasSize(3);

        assertTrue(sets.get(0)).containsExactlyInAnyOrder(one, three);
        assertTrue(sets.get(1)).containsExactlyInAnyOrder(two);
        assertTrue(sets.get(2)).containsExactlyInAnyOrder(four);
    }

    @Test
    public void canPhaseEvolving()
    {
        final SupportedAssembly one = assembly("A", "B", "C");
        final SupportedAssembly two = assembly("D", "E", "F");
        final SupportedAssembly three = assembly("A", "B", "C");
        final SupportedAssembly four = assembly("H");
        final SupportedAssembly five = assembly("A", "B", "H");

        final List<Set<SupportedAssembly>> sets = PrimaryPhasing.run(List.of(one, two, three, four, five));
        assertTrue(sets).hasSize(2);

        assertTrue(sets.get(0)).containsExactlyInAnyOrder(one, three, four, five);
        assertTrue(sets.get(1)).containsExactlyInAnyOrder(two);
    }

    @Test
    public void canPhaseSingle()
    {
        final SupportedAssembly one = assembly("A", "B", "C");
        final SupportedAssembly two = assembly("D", "E", "F");
        final SupportedAssembly three = assembly("A", "B", "C");
        final SupportedAssembly four = assembly("H");
        final SupportedAssembly five = assembly("A", "B", "H");
        final SupportedAssembly six = assembly("I", "D", "H");

        final List<Set<SupportedAssembly>> sets = PrimaryPhasing.run(List.of(one, two, three, four, five, six));
        assertTrue(sets).hasSize(1);

        assertTrue(sets.get(0)).containsExactlyInAnyOrder(one, two, three, four, five, six);
    }
    */
}