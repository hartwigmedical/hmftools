package com.hartwig.hmftools.lilac.coverage;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaComplexCoverageRanking
{
    private final int mMaxDistanceFromTopScore;
    private final List<HlaAllele> mCommon;
    private final List<HlaAllele> mRecovered;
    private final List<HlaAllele> mFavourites;

    public HlaComplexCoverageRanking(
            int maxDistanceFromTopScore, final List<HlaAllele> common, final List<HlaAllele> recovered,
            final List<HlaAllele> favourites)
    {
        mMaxDistanceFromTopScore = maxDistanceFromTopScore;
        mCommon = common;
        mRecovered = recovered;
        mFavourites = favourites;
    }

    public static List<HlaComplexCoverage> candidateRanking(final List<HlaComplexCoverage> complexes)
    {
        // TODO
        return Lists.newArrayList();
        
        /*
        int topScore;
        HlaComplexCoverage it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        boolean bl;
        Intrinsics.checkParameterIsNotNull(complexes, (String) "complexes");
        Collection collection = complexes;
        boolean bl2 = bl = !collection.isEmpty();
        if(!bl)
        {
            String string = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(string.toString());
        }
        Iterable iterable = $receiver$iv = (Iterable) complexes;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            HlaComplexCoverage
                    hlaComplexCoverage = (HlaComplexCoverage) item$iv$iv;
            Collection collection2 = destination$iv$iv;
            boolean bl3 = false;
            Integer n = it.getTotalCoverage();
            collection2.add(n);
        }
        Comparable comparable = CollectionsKt.max((Iterable) ((List) destination$iv$iv));
        if(comparable == null)
        {
            Intrinsics.throwNpe();
        }
        if((topScore = ((Number) ((Object) comparable)).intValue()) == 0)
        {
            return CollectionsKt.emptyList();
        }
        $receiver$iv = complexes;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaComplexCoverage) element$iv$iv;
            boolean bl4 = false;
            if(!(it.getTotalCoverage() >= topScore - mMaxDistanceFromTopScore))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        return CollectionsKt.sortedWith((Iterable) ((List) destination$iv$iv), (Comparator) new Comparator<HlaComplexCoverage>(this)
         */
    }

    private static int compare(HlaComplexCoverage o1, HlaComplexCoverage o2)
    {
        return 1;
        
        /*
        int o2RecoveredCount;
        int o2CommonCount;
        int favouritesCount = Intrinsics.compare((int) favourites(o1), (int) favourites(o2));
        if(favouritesCount != 0)
        {
            return -favouritesCount;
        }
        int wildcardCount = Intrinsics.compare((int) wildcardCount(o1), (int) wildcardCount(o2));
        if(wildcardCount != 0)
        {
            return wildcardCount;
        }
        int homozygousCompare = Intrinsics.compare((int) o1.homozygousAlleles(), (int) o2.homozygousAlleles());
        if(homozygousCompare != 0)
        {
            return -homozygousCompare;
        }
        int o1CommonCount = commonCount(o1);
        int commonCountCompare = Intrinsics.compare((int) o1CommonCount, (int) (o2CommonCount = commonCount(o2)));
        if(commonCountCompare != 0)
        {
            return -commonCountCompare;
        }
        int o1RecoveredCount = recoveredCount(o1);
        int recoveredCountCompare = Intrinsics.compare((int) o1RecoveredCount, (int) (o2RecoveredCount = recoveredCount(o2)));
        if(recoveredCountCompare != 0)
        {
            return recoveredCountCompare;
        }
        int n = 0;
        int n2 = o1.getAlleleCoverage().size();
        int n3 = o2.getAlleleCoverage().size();
        int n4 = Math.min(n2, n3);
        while(n < n4)
        {
            HlaAllele o2Allele;
            void i;
            HlaAllele o1Allele = o1.getAlleleCoverage().get((int) i).getAllele();
            int alleleCompare = o1Allele.compareTo(o2Allele = o2.getAlleleCoverage().get((int) i).getAllele());
            if(alleleCompare != 0)
            {
                return alleleCompare;
            }
            ++i;
        }
        throw (Throwable) new UnsupportedOperationException("Should not be able to make it to here");
        
         */
    }

    private static int favourites(final HlaComplexCoverage coverage)
    {
        return 1;
        
        /*
        HlaAllele it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) coverage.getAlleleCoverage();
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            HlaAlleleCoverage
                    hlaAlleleCoverage = (HlaAlleleCoverage) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaAllele hlaAllele = ((HlaAlleleCoverage) ((Object) it)).getAllele();
            collection.add(hlaAllele);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!mFavourites.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection = (List) destination$iv$iv;
        return collection.size();
        
         */
    }

    private static int wildcardCount(HlaComplexCoverage coverage)
    {
        return 1;
        // TODO
        
        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) coverage.getAlleleCoverage();
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            HlaAlleleCoverage it = (HlaAlleleCoverage) element$iv$iv;
            boolean bl = false;
            if(!(it.getWildCoverage() > 0.0))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection = (List) destination$iv$iv;
        return collection.size();
        
         */
    }

    private final int commonCount(HlaComplexCoverage coverage)
    {
        return 1;
        
        /*
        HlaAllele it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) coverage.getAlleleCoverage();
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            HlaAlleleCoverage
                    hlaAlleleCoverage = (HlaAlleleCoverage) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaAllele hlaAllele = ((HlaAlleleCoverage) ((Object) it)).getAllele();
            collection.add(hlaAllele);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!mCommon.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection = (List) destination$iv$iv;
        return collection.size();
        
         */
    }

    private final int recoveredCount(HlaComplexCoverage coverage)
    {
        return 1;
        
        /*
        HlaAllele it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) coverage.getAlleleCoverage();
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            HlaAlleleCoverage
                    hlaAlleleCoverage = (HlaAlleleCoverage) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaAllele hlaAllele = ((HlaAlleleCoverage) ((Object) it)).getAllele();
            collection.add(hlaAllele);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (HlaAllele) element$iv$iv;
            boolean bl = false;
            if(!mRecovered.contains(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        Collection collection = (List) destination$iv$iv;
        return collection.size();
        
         */
    }
}
