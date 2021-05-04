package com.hartwig.hmftools.lilac.coverage;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

public class HlaComplex
{
    private final List<HlaAllele> Alleles;

    public HlaComplex(final List<HlaAllele> alleles)
    {
        Alleles = alleles;
    }

    public static List<HlaComplex> complexes(final LilacConfig config, final List<FragmentAlleles> referenceFragmentAlleles,
            final List<HlaAllele> candidateAlleles, final List<HlaAllele> recoveredAlleles)
    {
        // TODO
        return Lists.newArrayList();
//        List<HlaComplex> list;
//        Collection<HlaAllele> collection;
//        void $receiver$iv$iv;
//        void $receiver$iv;
//        Iterable $receiver$iv$iv2;
//        Iterable $receiver$iv2;
//        void $receiver$iv$iv3;
//        Iterable $receiver$iv3;
//        Intrinsics.checkParameterIsNotNull((Object) config, (String) "config");
//        Intrinsics.checkParameterIsNotNull(referenceFragmentAlleles, (String) "referenceFragmentAlleles");
//        Intrinsics.checkParameterIsNotNull(candidateAlleles, (String) "candidateAlleles");
//        Intrinsics.checkParameterIsNotNull(recoveredAlleles, (String) "recoveredAlleles");
//        getLogger().info("Identifying uniquely identifiable groups and proteins [total,unique,shared,wide]");
//        HlaComplexCoverage groupCoverage =
//                HlaComplexCoverageFactory.Companion.groupCoverage(referenceFragmentAlleles, (Collection<HlaAllele>) candidateAlleles);
//        List<HlaAlleleCoverage> confirmedGroups = confirmUnique(groupCoverage, config);
//        Iterable iterable = $receiver$iv3 = (Iterable) groupCoverage.getAlleleCoverage();
//        Collection destination$iv$iv = new ArrayList();
//        for(Object element$iv$iv : $receiver$iv$iv3)
//        {
//            HlaAlleleCoverage it = (HlaAlleleCoverage) element$iv$iv;
//            boolean bl = false;
//            if(!(it.getUniqueCoverage() > 0 && !confirmedGroups.contains(it)))
//            {
//                continue;
//            }
//            destination$iv$iv.add(element$iv$iv);
//        }
//        List discardedGroups = CollectionsKt.sortedDescending((Iterable) ((List) destination$iv$iv));
//        $receiver$iv3 = confirmedGroups;
//        if(!$receiver$iv3.isEmpty())
//        {
//            getLogger()
//                    .info("    confirmed " + confirmedGroups.size() + " unique groups: "
//                            + CollectionsKt.joinToString$default((Iterable) confirmedGroups, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
//        }
//        else
//        {
//            getLogger().info("    confirmed 0 unique groups");
//        }
//        $receiver$iv3 = discardedGroups;
//        if(!$receiver$iv3.isEmpty())
//        {
//            getLogger()
//                    .info("    found " + discardedGroups.size() + " insufficiently unique groups: "
//                            + CollectionsKt.joinToString$default((Iterable) discardedGroups, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
//        }
//        List<HlaAllele> confirmedGroupAlleles = alleles(confirmedGroups);
//        List<HlaAllele> candidatesAfterConfirmedGroups = filterWithConfirmedGroups(candidateAlleles, confirmedGroupAlleles);
//        HlaComplexCoverage proteinCoverage =
//                HlaComplexCoverageFactory.Companion.proteinCoverage(referenceFragmentAlleles, (Collection<HlaAllele>) candidatesAfterConfirmedGroups);
//        List<HlaAlleleCoverage> confirmedProtein = confirmUnique(proteinCoverage, config);
//        Iterable bl = $receiver$iv2 = (Iterable) proteinCoverage.getAlleleCoverage();
//        Collection destination$iv$iv2 = new ArrayList();
//        for(Object element$iv$iv : $receiver$iv$iv2)
//        {
//            HlaAlleleCoverage it = (HlaAlleleCoverage) element$iv$iv;
//            boolean bl2 = false;
//            if(!(it.getUniqueCoverage() > 0 && !confirmedProtein.contains(it)))
//            {
//                continue;
//            }
//            destination$iv$iv2.add(element$iv$iv);
//        }
//        List discardedProtein = CollectionsKt.sortedDescending((Iterable) ((List) destination$iv$iv2));
//        $receiver$iv2 = confirmedProtein;
//        if(!$receiver$iv2.isEmpty())
//        {
//            getLogger()
//                    .info("    confirmed " + confirmedProtein.size() + " unique proteins: "
//                            + CollectionsKt.joinToString$default((Iterable) confirmedProtein, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
//        }
//        else
//        {
//            getLogger().info("    confirmed 0 unique proteins");
//        }
//        $receiver$iv2 = discardedProtein;
//        if(!$receiver$iv2.isEmpty())
//        {
//            getLogger()
//                    .info("    found " + discardedProtein.size() + " insufficiently unique proteins: "
//                            + CollectionsKt.joinToString$default((Iterable) discardedProtein, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
//        }
//        $receiver$iv$iv2 = confirmedProtein;
//        List<HlaAllele> list2 = candidatesAfterConfirmedGroups;
//        Companion companion = this;
//        destination$iv$iv2 = $receiver$iv;
//        Collection destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
//        for(Object item$iv$iv : $receiver$iv$iv)
//        {
//            void it;
//            HlaAlleleCoverage bl2 = (HlaAlleleCoverage) item$iv$iv;
//            collection = destination$iv$iv3;
//            boolean bl3 = false;
//            HlaAllele hlaAllele = it.getAllele();
//            collection.add(hlaAllele);
//        }
//        collection = (List) destination$iv$iv3;
//        List<HlaAllele> candidatesAfterConfirmedProteins = companion.filterWithConfirmedProteins(list2, (List<HlaAllele>) collection);
//        List<HlaComplex> aOnlyComplexes =
//                gene("A", alleles(confirmedGroups), alleles(confirmedProtein), candidatesAfterConfirmedProteins);
//        List<HlaComplex> bOnlyComplexes =
//                gene("B", alleles(confirmedGroups), alleles(confirmedProtein), candidatesAfterConfirmedProteins);
//        List<HlaComplex> cOnlyComplexes =
//                gene("C", alleles(confirmedGroups), alleles(confirmedProtein), candidatesAfterConfirmedProteins);
//        List<HlaComplex> complexes = null;
//        long simpleComplexCount = (long) aOnlyComplexes.size() * (long) bOnlyComplexes.size() * (long) cOnlyComplexes.size();
//        if(simpleComplexCount > (long) 100000 || simpleComplexCount < 0L)
//        {
//            getLogger().info("Candidate permutations exceeds maximum complexity");
//            HlaComplexCoverageFactory groupRankedCoverageFactory = new HlaComplexCoverageFactory(config);
//            List<HlaAllele> aTopCandidates =
//                    groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, aOnlyComplexes);
//            List<HlaAllele> bTopCandidates =
//                    groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, bOnlyComplexes);
//            List<HlaAllele> cTopCandidates =
//                    groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, cOnlyComplexes);
//            List topCandidates =
//                    CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) aTopCandidates, (Iterable) bTopCandidates), (Iterable) cTopCandidates);
//            Set rejected = CollectionsKt.subtract((Iterable) candidatesAfterConfirmedProteins, (Iterable) topCandidates);
//            getLogger()
//                    .info("    discarding " + rejected.size() + " unlikely candidates: "
//                            + CollectionsKt.joinToString$default((Iterable) rejected, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
//            list = complexes(alleles(confirmedGroups), alleles(confirmedProtein), topCandidates);
//        }
//        else
//        {
//            list = complexes(alleles(confirmedGroups), alleles(confirmedProtein), candidatesAfterConfirmedProteins);
//        }
//        complexes = list;
//        return complexes;
    }

    private final List<HlaComplex> complexes(List<HlaAllele> confirmedGroups, List<HlaAllele> confirmedProteins,
            List<HlaAllele> candidates)
    {
        List<HlaComplex> a = gene("A", confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> b = gene("B", confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> c = gene("C", confirmedGroups, confirmedProteins, candidates);
        return combineComplexes(combineComplexes(a, b), c);
    }

    /*
     * WARNING - void declaration
     */
    public static List<HlaComplex> gene(final String gene, final List<HlaAllele> unfilteredGroups,
            final List<HlaAllele> unfilteredProteins, final List<HlaAllele> unfilteredCandidates)
    {
        // TODO
        return Lists.newArrayList();

    }

    public static List<HlaComplex> combineComplexes(final List<HlaComplex> first, final List<HlaComplex> second)
    {
        // TODO
        return Lists.newArrayList();
    }

    private static List<HlaComplex> combineAlleles(List<HlaAllele> first, List<HlaAllele> second)
    {
        // TODO
        return Lists.newArrayList();
        
    }

    /*
    private static <T> List<List<T>> cartesianProduct(Collection<? extends T> first, Collection<? extends T> second)
    {
        List result = new ArrayList();
        for(T i : first)
        {
            for(T j : second)
            {
                if(!(Intrinsics.areEqual(i, j) ^ true))
                {
                    continue;
                }
                result.add(CollectionsKt.mutableListOf((Object[]) new Object[] { i, j }));
            }
        }
        return CollectionsKt.distinct((Iterable) result);
    }
    
     */

    private static List<HlaAllele> filterWithConfirmedProteins(List<HlaAllele> $receiver, List<HlaAllele> confirmedGroups)
    {
        // TODO
        return Lists.newArrayList();
    }

    private static List<HlaAllele> filterWithConfirmedGroups(List<HlaAllele> $receiver, List<HlaAllele> confirmedGroups)
    {
        // TODO
        return Lists.newArrayList();
    }

    private final List<HlaAllele> filterWithConfirmed(List<HlaAllele> $receiver, List<HlaAllele> confirmed) //            Function1<? super HlaAllele, HlaAllele> transform)
    {
        // TODO
        return Lists.newArrayList();
    }

    private static List<HlaAlleleCoverage> confirmUnique(final HlaComplexCoverage coverage, final LilacConfig config)
    {
        // TODO
        return Lists.newArrayList();
    }

    private static List<HlaAllele> alleles(List<HlaAlleleCoverage> $receiver)
    {
        // TODO
        return Lists.newArrayList();
    }
}
