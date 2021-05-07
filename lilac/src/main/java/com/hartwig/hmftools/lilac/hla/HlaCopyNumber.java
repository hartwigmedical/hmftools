package com.hartwig.hmftools.lilac.hla;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public final class HlaCopyNumber
{
    public final HlaAllele Allele;
    public final double CopyNumber;

    public HlaCopyNumber(final HlaAllele allele, double copyNumber)
    {
        Allele = allele;
        CopyNumber = copyNumber;
    }

    public static List<HlaCopyNumber> alleleCopyNumber(final List<HlaAllele> winners)
    {
        return winners.stream().map(x -> new HlaCopyNumber(x, 0)).collect(Collectors.toList());
    }

    public static List<HlaCopyNumber> alleleCopyNumber(
            final List<HlaAllele> winners, final String geneCopyNumberFile, final HlaComplexCoverage tumorCoverage)
    {
        /*
                    if (geneCopyNumberFile.isEmpty() || tumorCoverage.alleleCoverage.isEmpty()) {
                return winners.map { HlaCopyNumber(it, 0.0) }
            }

            val coverage = tumorCoverage.alleleCoverage
            require(coverage.size == 6)
            val geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberFile).filter { it.gene() in genes}.sortedBy { it.gene() }
            val aCopyNumber = alleleCopyNumber(geneCopyNumbers[0], coverage.filter { it.allele.gene == "A" })
            val bCopyNumber = alleleCopyNumber(geneCopyNumbers[1], coverage.filter { it.allele.gene == "B" })
            val cCopyNumber = alleleCopyNumber(geneCopyNumbers[2], coverage.filter { it.allele.gene == "C" })

            return (aCopyNumber + bCopyNumber + cCopyNumber).sortedBy { it.allele }

         */

        if(geneCopyNumberFile.isEmpty() || tumorCoverage.getAlleleCoverage().isEmpty())
            return alleleCopyNumber(winners);


        return Lists.newArrayList();

        /*
        CharSequence charSequence = geneCopyNumberFile;
        if(charSequence.length() == 0 || tumorCoverage.getAlleleCoverage().isEmpty())
        {
            void $receiver$iv$iv52;
            Iterable $receiver$iv5;
            Iterable iterable = $receiver$iv5 = (Iterable) winners;
            Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv5, (int) 10));
            for(Object item$iv$iv : $receiver$iv$iv52)
            {
                void it;
                HlaAllele hlaAllele = (HlaAllele) item$iv$iv;
                Collection collection = destination$iv$iv;
                boolean bl = false;
                HlaCopyNumber hlaCopyNumber = new HlaCopyNumber((HlaAllele) it, 0.0);
                collection.add(hlaCopyNumber);
            }
            return (List) destination$iv$iv;
        }
        List<com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage> coverage = tumorCoverage.getAlleleCoverage();
        boolean bl = $receiver$iv$iv52 = coverage.size() == 6;
        if(!$receiver$iv$iv52)
        {
            String destination$iv$iv = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(destination$iv$iv.toString());
        }
        List list = GeneCopyNumberFile.read((String) geneCopyNumberFile);
        Intrinsics.checkExpressionValueIsNotNull((Object) list, (String) "GeneCopyNumberFile.read(geneCopyNumberFile)");
        Iterable iterable = $receiver$iv4 = (Iterable) list;
        Object destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv4)
        {
            GeneCopyNumber it = (GeneCopyNumber) element$iv$iv;
            boolean bl2 = false;
            if(!genes.contains(it.gene()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv4 = (List) destination$iv$iv;
        $receiver$iv$iv4 = $receiver$iv4;
        destination$iv$iv = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                GeneCopyNumber it = (GeneCopyNumber) a;
                boolean bl = false;
                Comparable comparable = (Comparable) ((Object) it.gene());
                it = (GeneCopyNumber) b;
                Comparable comparable2 = comparable;
                bl = false;
                String string = it.gene();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) ((Comparable) ((Object) string)));
            }
        };
        List geneCopyNumbers = CollectionsKt.sortedWith((Iterable) $receiver$iv$iv4, (Comparator) destination$iv$iv);
        Object e = geneCopyNumbers.get(0);
        Intrinsics.checkExpressionValueIsNotNull(e, (String) "geneCopyNumbers[0]");
        $receiver$iv$iv4 = coverage;
        GeneCopyNumber geneCopyNumber = (GeneCopyNumber) e;
        Companion companion = this;
        destination$iv$iv = $receiver$iv3;
        Collection destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv3)
        {
            com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage it = (com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage) element$iv$iv;
            boolean bl3 = false;
            if(!Intrinsics.areEqual((Object) it.getAllele().getGene(), (Object) "A"))
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        List list2 = (List) destination$iv$iv2;
        List<HlaCopyNumber> aCopyNumber = companion.alleleCopyNumber(geneCopyNumber, list2);
        Object e2 = geneCopyNumbers.get(1);
        Intrinsics.checkExpressionValueIsNotNull(e2, (String) "geneCopyNumbers[1]");
        $receiver$iv$iv3 = coverage;
        geneCopyNumber = (GeneCopyNumber) e2;
        companion = this;
        destination$iv$iv2 = $receiver$iv2;
        Collection destination$iv$iv3 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage it = (com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage) element$iv$iv;
            boolean bl4 = false;
            if(!Intrinsics.areEqual((Object) it.getAllele().getGene(), (Object) "B"))
            {
                continue;
            }
            destination$iv$iv3.add(element$iv$iv);
        }
        list2 = (List) destination$iv$iv3;
        List<HlaCopyNumber> bCopyNumber = companion.alleleCopyNumber(geneCopyNumber, list2);
        Object e3 = geneCopyNumbers.get(2);
        Intrinsics.checkExpressionValueIsNotNull(e3, (String) "geneCopyNumbers[2]");
        $receiver$iv$iv2 = coverage;
        geneCopyNumber = (GeneCopyNumber) e3;
        companion = this;
        destination$iv$iv3 = $receiver$iv;
        Collection destination$iv$iv4 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage it = (com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage) element$iv$iv;
            boolean bl5 = false;
            if(!Intrinsics.areEqual((Object) it.getAllele().getGene(), (Object) "C"))
            {
                continue;
            }
            destination$iv$iv4.add(element$iv$iv);
        }
        list2 = (List) destination$iv$iv4;
        List<HlaCopyNumber> cCopyNumber = companion.alleleCopyNumber(geneCopyNumber, list2);
        Iterable iterable2 = $receiver$iv =
                (Iterable) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) aCopyNumber, (Iterable) bCopyNumber), (Iterable) cCopyNumber);
        Comparator comparator = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                HlaCopyNumber it = (HlaCopyNumber) a;
                boolean bl = false;
                Comparable comparable = it.getAllele();
                it = (HlaCopyNumber) b;
                Comparable comparable2 = comparable;
                bl = false;
                HlaAllele hlaAllele = it.getAllele();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) hlaAllele);
            }
        };
        return CollectionsKt.sortedWith((Iterable) iterable2, (Comparator) comparator);

         */
    }

    private static List<HlaCopyNumber> alleleCopyNumber(GeneCopyNumber geneCopyNumber, List<HlaAlleleCoverage> alleleCoverage)
    {
        return Lists.newArrayList();
        // TODO

        /*
        boolean bl;
        boolean bl2 = bl = alleleCoverage.size() == 2;
        if(!bl)
        {
            String string = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(string.toString());
        }
        double minor = geneCopyNumber.minMinorAlleleCopyNumber();
        double major = geneCopyNumber.minCopyNumber() - minor;
        return alleleCoverage.get(0).getTotalCoverage() >= alleleCoverage.get(1).getTotalCoverage()
                ? CollectionsKt.listOf((Object[]) new HlaCopyNumber[] { new HlaCopyNumber(alleleCoverage.get(0).getAllele(), major),
                new HlaCopyNumber(alleleCoverage.get(1).getAllele(), minor) })
                : CollectionsKt.listOf((Object[]) new HlaCopyNumber[] { new HlaCopyNumber(alleleCoverage.get(0).getAllele(), minor),
                        new HlaCopyNumber(alleleCoverage.get(1).getAllele(), major) });

         */
    }

}
