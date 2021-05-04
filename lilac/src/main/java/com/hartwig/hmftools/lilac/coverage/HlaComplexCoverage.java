package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.round;

import java.io.File;
import java.util.List;

public final class HlaComplexCoverage implements Comparable<HlaComplexCoverage>
{
    public final int TotalCoverage;
    public final int UniqueCoverage;
    public final int SharedCoverage;
    public final int WildCoverage;
    
    public final List<HlaAlleleCoverage> AlleleCoverage;

    public HlaComplexCoverage(int uniqueCoverage, int sharedCoverage, int wildCoverage, final List<HlaAlleleCoverage> alleleCoverage)
    {
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        AlleleCoverage = alleleCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;
    }

    public final HlaComplexCoverage expandToSixAlleles()
    {
        return create(HlaAlleleCoverage.expand(AlleleCoverage));
    }

    public final int homozygousAlleles()
    {
        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) AlleleCoverage;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            HlaAlleleCoverage
                    hlaAlleleCoverage = (HlaAlleleCoverage) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaAllele hlaAllele = it.getAllele();
            collection.add(hlaAllele);
        }
        int distinct = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv)).size();
        return 6 - distinct;

         */

        return 0;
    }

    @Override
    public int compareTo(final HlaComplexCoverage other)
    {
        /*
        int totalCoverageCompare = Intrinsics.compare((int) TotalCoverage, (int) other.TotalCoverage);
        if(totalCoverageCompare != 0)
        {
            return totalCoverageCompare;
        }
        int wildCoverageCompare = Intrinsics.compare((int) WildCoverage, (int) other.WildCoverage);
        if(wildCoverageCompare != 0)
        {
            return -wildCoverageCompare;
        }
        int sharedCoverageCompare = Intrinsics.compare((int) SharedCoverage, (int) other.SharedCoverage);
        if(sharedCoverageCompare != 0)
        {
            return sharedCoverageCompare;
        }
        return Intrinsics.compare((int) UniqueCoverage, (int) other.UniqueCoverage);

         */

        // TODO
        return 1;
    }

    public String toString()
    {
        //         val types = alleleCoverage.map { it.allele }.toSet()
        //        return "${totalCoverage}\t$uniqueCoverage\t${sharedCoverage}\t${wildCoverage}\t${types.size}\t${alleleCoverage.joinToString("\t")}"
        return "" ;
        // TODO
    }

    public final List<HlaAlleleCoverage> getAlleleCoverage()
    {
        return AlleleCoverage;
    }

    public static HlaComplexCoverage create(final List<HlaAlleleCoverage> alleles)
    {
        int unique = 0;
        double shared = 0.0;
        double wild = 0.0;
        for(HlaAlleleCoverage coverage : alleles)
        {
            unique += coverage.UniqueCoverage;
            shared += coverage.SharedCoverage;
            wild += coverage.WildCoverage;
        }

        // TODO
        // return new HlaComplexCoverage(unique, (int)round(shared), (int)round(wild), alleles.sortedBy { it.allele })
        return null;

    }

    public final String header()
    {
        return "totalCoverage\tuniqueCoverage\tsharedCoverage\twildCoverage\ttypes\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6";
    }

    public final void writeToFile(final List<HlaComplexCoverage> coverage, final String fileName)
    {
        File file = new File(fileName);

        /*
        FilesKt.writeText$default((File) file, (String) (header() + "\n"), null, (int) 2, null);
        for(HlaComplexCoverage coverage : $receiver)
        {
            FilesKt.appendText$default((File) file, (String) (coverage.toString() + "\n"), null, (int) 2, null);
        }
        */
    }
}
