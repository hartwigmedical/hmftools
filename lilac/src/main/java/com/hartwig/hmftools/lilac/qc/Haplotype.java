package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacUtils.listMax;
import static com.hartwig.hmftools.lilac.LilacUtils.listMin;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.lilac.SequenceCount;
import org.apache.commons.math3.util.Pair;

public class Haplotype
{
    public final int StartLocus;
    public final int EndLocus;
    public final int SupportingFragments;
    public final String Haplotype;

    public Haplotype(final int startLocus, final int endLocus, final int supportingFragments, final String haplotype)
    {
        StartLocus = startLocus;
        EndLocus = endLocus;
        SupportingFragments = supportingFragments;
        Haplotype = haplotype;
    }

    public static Haplotype fromString(final String line)
    {
        String[] items = line.split(DELIM, -1);
        int startLocus = Integer.parseInt(items[0]);
        String haplotype = items[1];
        int endLocus = startLocus + haplotype.length() - 1;
        return new Haplotype(startLocus, endLocus, 0, haplotype);
    }

    public static Haplotype create(final List<Integer> aminoAcidIndices, final Pair<String,Integer> evidence, final SequenceCount aminoAcidCount)
    {
        int startLocus = listMin(aminoAcidIndices);
        int endLocus = listMax(aminoAcidIndices);
        String sparseHaplotype = evidence.getFirst();
        String complexHaplotype = "";

        for(int locus = startLocus; locus <= endLocus; ++locus)
        {
            if(aminoAcidIndices.contains(locus))
                complexHaplotype += sparseHaplotype.charAt(aminoAcidIndices.indexOf(locus));
            else
                complexHaplotype += aminoAcidCount.getMinCountSequences(locus).get(0);
        }

        return new Haplotype(startLocus, endLocus, evidence.getSecond(), complexHaplotype);
    }

    public boolean contains(final Haplotype unmatched)
    {
        if (unmatched.StartLocus >= StartLocus && unmatched.EndLocus <= EndLocus) {
            int start = max(unmatched.StartLocus, StartLocus);
            int end = min(unmatched.EndLocus, EndLocus);

            // CHECK off by
            for(int locus = start; locus < end; ++locus) {
                if(charAt(locus) != unmatched.charAt(locus))
                    return false;
            }

            return true;
        }

        return false;
    }

    private char charAt(int locus)
    {
        int index = locus - StartLocus;
        return Haplotype.charAt(index);
    }

    public String toString()
    {
        return String.format("startLocus=%d, endLocus=%d, supportingFragments=%d, haplotype=%s",
                StartLocus, EndLocus, SupportingFragments, Haplotype);
    }

    public static class HaplotypeFragmentSorter implements Comparator<Haplotype>
    {
        public int compare(final Haplotype first, final Haplotype second)
        {
            int compare = first.SupportingFragments - second.SupportingFragments;
            if(compare != 0)
                return compare > 0 ? 1 : -1;

            return 0;
        }
    }

    public static class HaplotypeStartLocusSorter implements Comparator<Haplotype>
    {
        public int compare(final Haplotype first, final Haplotype second)
        {
            int compare = first.StartLocus - second.StartLocus;
            if(compare != 0)
                return compare < 0 ? 1 : -1; // lower start first

            return 0;
        }

    }

}
