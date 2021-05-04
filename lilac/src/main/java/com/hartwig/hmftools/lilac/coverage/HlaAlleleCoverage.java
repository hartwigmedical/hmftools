package com.hartwig.hmftools.lilac.coverage;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

import org.jetbrains.annotations.NotNull;

public class HlaAlleleCoverage implements Comparable<HlaAlleleCoverage>
{
    public final HlaAllele Allele;
    public final int UniqueCoverage;
    public final double SharedCoverage;
    public final double WildCoverage;
    public final double TotalCoverage;

    public HlaAlleleCoverage(final HlaAllele allele, int uniqueCoverage, double sharedCoverage, double wildCoverage)
    {
        Allele = allele;
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;
    }

    public int compareTo(final HlaAlleleCoverage other)
    {
        // TODO

        /*
        int uniqueCompare = Intrinsics.compare((int) UniqueCoverage, (int) other.UniqueCoverage);
        if(uniqueCompare != 0)
        {
            return uniqueCompare;
        }
        return Double.compare(SharedCoverage, other.SharedCoverage);

         */
        return 1;
    }

    @NotNull
    public String toString()
    {
        // "totalCoverage\tuniqueCoverage\tsharedCoverage\twildCoverage\ttypes\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6"
        return "";
    }

    public static List<HlaAlleleCoverage> expand(@NotNull List<HlaAlleleCoverage> $receiver)
    {
        return Lists.newArrayList();
    }

    public static List<HlaAlleleCoverage> proteinCoverage(@NotNull List<FragmentAlleles> fragmentSequences)
    {
        // return create(fragmentSequences) { it }
        // TODO
        return Lists.newArrayList();
    }

    public static List<HlaAlleleCoverage> groupCoverage(@NotNull List<FragmentAlleles> fragmentSequences)
    {
        // TODO
        // return create(fragmentSequences) { it.asAlleleGroup() }
        return Lists.newArrayList();
    }

    public static List<HlaAlleleCoverage> create(@NotNull List<FragmentAlleles> fragmentSequences) // pass in functor? final Function1<? super HlaAllele, HlaAllele> type
    {
        // TODO
        // see orig
        return Lists.newArrayList();
    }

    private static List<HlaAlleleCoverage> splitSingle(List<HlaAlleleCoverage> $receiver)
    {
        // TODO
        // see orig
        return Lists.newArrayList();
    }
}
