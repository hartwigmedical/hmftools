package com.hartwig.hmftools.svassembly.processor;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.svassembly.models.SupportedAssembly;

import org.junit.Test;

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

    @Test
    public void canPhaseSimple()
    {
        final SupportedAssembly one = assembly("A", "B", "C");
        final SupportedAssembly two = assembly("D", "E", "F");
        final SupportedAssembly three = assembly("A", "B", "C");
        final SupportedAssembly four = assembly("H");

        final List<Set<SupportedAssembly>> sets = PrimaryPhasing.run(List.of(one, two, three, four));
        assertThat(sets).hasSize(3);

        assertThat(sets.get(0)).containsExactlyInAnyOrder(one, three);
        assertThat(sets.get(1)).containsExactlyInAnyOrder(two);
        assertThat(sets.get(2)).containsExactlyInAnyOrder(four);
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
        assertThat(sets).hasSize(2);

        assertThat(sets.get(0)).containsExactlyInAnyOrder(one, three, four, five);
        assertThat(sets.get(1)).containsExactlyInAnyOrder(two);
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
        assertThat(sets).hasSize(1);

        assertThat(sets.get(0)).containsExactlyInAnyOrder(one, two, three, four, five, six);
    }
}