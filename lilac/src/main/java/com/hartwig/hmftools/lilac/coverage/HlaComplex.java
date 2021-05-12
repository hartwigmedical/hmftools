package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.groupCoverage;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.proteinCoverage;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.contains;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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

    public List<HlaAllele> getAlleles() { return Alleles; }

    public static List<HlaComplex> complexes(
            final LilacConfig config, final List<FragmentAlleles> referenceFragmentAlleles,
            final List<HlaAllele> candidateAlleles) // final List<HlaAllele> recoveredAlleles was uused
    {
        LL_LOGGER.info("Identifying uniquely identifiable groups and proteins [total,unique,shared,wild]");

        HlaComplexCoverage groupCoverage = groupCoverage(referenceFragmentAlleles, candidateAlleles);

        List<HlaAlleleCoverage> confirmedGroups = groupCoverage.confirmUnique(config);

        List<HlaAlleleCoverage> discardedGroups = groupCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !confirmedGroups.contains(x)).collect(Collectors.toList());

        Collections.sort(discardedGroups, Collections.reverseOrder());

        if(!confirmedGroups.isEmpty())
        {
            LL_LOGGER.info("  confirmed {} unique groups: {}",  confirmedGroups.size(), HlaAlleleCoverage.toString(confirmedGroups));
        }
        else
        {
            LL_LOGGER.info("  confirmed 0 unique groups");
        }

        if (!discardedGroups.isEmpty())
        {
            LL_LOGGER.info("  found {} insufficiently unique groups: {}",
                    discardedGroups.size(), HlaAlleleCoverage.toString(discardedGroups));
        }

        List<HlaAllele> confirmedGroupAlleles = alleles(confirmedGroups);
        List<HlaAllele> candidatesAfterConfirmedGroups = filterWithConfirmedGroups(candidateAlleles, confirmedGroupAlleles);
        HlaComplexCoverage proteinCoverage = proteinCoverage(referenceFragmentAlleles, candidatesAfterConfirmedGroups);
        List<HlaAlleleCoverage> confirmedProtein = proteinCoverage.confirmUnique(config);
        List<HlaAlleleCoverage> discardedProtein = proteinCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !confirmedProtein.contains(x)).collect(Collectors.toList());
        Collections.sort(discardedProtein, Collections.reverseOrder());

        if(!confirmedProtein.isEmpty())
        {
            LL_LOGGER.info("  confirmed {} unique groups: {}", confirmedProtein.size(), HlaAlleleCoverage.toString(confirmedProtein));
        }
        else
        {
            LL_LOGGER.info("  confirmed 0 unique groups");
        }

        if (!discardedProtein.isEmpty())
        {
            LL_LOGGER.info("  found {} insufficiently unique groups: {}", discardedProtein.size(), HlaAlleleCoverage.toString(discardedProtein));
        }

        List<HlaAllele> confirmedProteinAlleles = alleles(confirmedProtein);

        List<HlaAllele> candidatesAfterConfirmedProteins = filterWithConfirmedProteins(candidatesAfterConfirmedGroups, confirmedProteinAlleles);

        List<HlaComplex> aOnlyComplexes = gene("A", confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);
        List<HlaComplex> bOnlyComplexes = gene("B", confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);
        List<HlaComplex> cOnlyComplexes = gene("C", confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);

        List<HlaComplex> complexes;
        long simpleComplexCount = aOnlyComplexes.size() * bOnlyComplexes.size() * cOnlyComplexes.size();

        if (simpleComplexCount > 100_000 || simpleComplexCount < 0)
        {
            LL_LOGGER.info("Candidate permutations exceeds maximum complexity");

            HlaComplexCoverageFactory groupRankedCoverageFactory = new HlaComplexCoverageFactory(config);
            List<HlaAllele> aTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, aOnlyComplexes);
            List<HlaAllele> bTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, bOnlyComplexes);
            List<HlaAllele> cTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, cOnlyComplexes);
            List<HlaAllele> topCandidates = Lists.newArrayList();
            topCandidates.addAll(aTopCandidates);
            topCandidates.addAll(bTopCandidates);
            topCandidates.addAll(cTopCandidates);

            List<HlaAllele> rejected = candidatesAfterConfirmedProteins.stream()
                    .filter(x -> !contains(topCandidates, x)).collect(Collectors.toList());

            LL_LOGGER.info("  discarding {} unlikely candidates: ", rejected.size(), HlaAllele.toString(rejected));
            complexes = buildComplexes(confirmedGroupAlleles, confirmedProteinAlleles, topCandidates);
        }
        else
        {
            complexes = buildComplexes(confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);
        }

        return complexes;
    }

    private static List<HlaComplex> buildComplexes(
            final List<HlaAllele> confirmedGroups, final List<HlaAllele> confirmedProteins, final List<HlaAllele> candidates)
    {
        List<HlaComplex> a = gene("A", confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> b = gene("B", confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> c = gene("C", confirmedGroups, confirmedProteins, candidates);
        return combineComplexes(combineComplexes(a, b), c);
    }

    public static List<HlaComplex> gene(
            final String gene, final List<HlaAllele> unfilteredGroups,
            final List<HlaAllele> unfilteredProteins, final List<HlaAllele> unfilteredCandidates)
    {
        List<HlaAllele> confirmedGroups = HlaAllele.takeN(unfilteredGroups.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList()), 2);
        List<HlaAllele> confirmedProteins = HlaAllele.takeN(unfilteredProteins.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList()), 2);
        List<HlaAllele> candidates = unfilteredCandidates.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList());

        if (confirmedProteins.size() == 2)
            return Lists.newArrayList(new HlaComplex(confirmedProteins));

        if (confirmedProteins.size() == 1)
        {
            List<HlaAllele> confirmedProteinGroups = confirmedProteins.stream().map(x -> x.asAlleleGroup()).collect(Collectors.toList());
            List<HlaAllele> remainingGroups = confirmedGroups.stream().filter(x -> !contains(confirmedProteinGroups, x)).collect(Collectors.toList());

            List<HlaAllele> first = confirmedProteins;

            List<HlaAllele> second = remainingGroups.isEmpty() ?
                    candidates.stream().filter(x -> !x.matches(confirmedProteins.get(0))).collect(Collectors.toList()) :
                    candidates.stream().filter(x -> contains(remainingGroups, x.asAlleleGroup())).collect(Collectors.toList());

            List<HlaComplex> complexes = combineAlleles(first, second);
            if(!remainingGroups.isEmpty())
                return complexes;

            // CHECK what does + do, combine the alleles?
            complexes.add(new HlaComplex(first));
            return complexes;

            // return if (remainingGroups.isEmpty()) combineAlleles(first, second) + HlaComplex(first) else combineAlleles(first, second)
        }

        if (confirmedGroups.size() == 2)
        {
            List<HlaAllele> first = candidates.stream().filter(x -> x.asAlleleGroup().matches(confirmedGroups.get(0))).collect(Collectors.toList());
            List<HlaAllele> second = candidates.stream().filter(x -> x.asAlleleGroup().matches(confirmedGroups.get(1))).collect(Collectors.toList());
            return combineAlleles(first, second);
        }

        if (confirmedGroups.size() == 1)
        {
            List<HlaAllele> first = candidates.stream().filter(x -> x.asAlleleGroup().matches(confirmedGroups.get(0))).collect(Collectors.toList());
            List<HlaAllele> second = candidates;

            List<HlaComplex> complexes = first.stream().map(x -> new HlaComplex(Lists.newArrayList(x))).collect(Collectors.toList());
            complexes.addAll(combineAlleles(first, second));
            return complexes;
        }

        List<HlaComplex> complexes = candidates.stream().map(x -> new HlaComplex(Lists.newArrayList(x))).collect(Collectors.toList());
        complexes.addAll(combineAlleles(candidates, candidates));
        return complexes;
    }

    public static List<HlaComplex> combineComplexes(final List<HlaComplex> first, final List<HlaComplex> second)
    {
        List<List<HlaComplex>> intermediatePairs = cartesianComplexProduct(first, second);
        // return intermediatePairs.stream().map(x -> )

        // TODO, check how combined
        // input is a list of 1 or 2 alleles, and then makes a new list of the 2 HlaComplex pairs (not combing them)`
        return Lists.newArrayList();
        //return intermediate.map { it.flatMap { it.alleles } }.map { HlaComplex(it) }
    }

    private static List<List<HlaComplex>> cartesianComplexProduct(final List<HlaComplex> first, final List<HlaComplex> second)
    {
        List<List<HlaComplex>> results = Lists.newArrayList();

        for(HlaComplex i : first)
        {
            for(HlaComplex j : second)
            {
                if(i != j)
                {
                    List<HlaComplex> pairing = Lists.newArrayList(i, j);

                    if(results.stream().anyMatch(x -> x.get(0) == pairing.get(0) && x.get(1) == pairing.get(1)))
                        continue;

                    results.add(pairing);
                }
            }
        }

        return results;
    }



    private static List<HlaComplex> combineAlleles(final List<HlaAllele> first, final List<HlaAllele> second)
    {
        List<List<HlaAllele>> allelePairs = cartesianAlleleProduct(first, second);
        return allelePairs.stream().map(x -> new HlaComplex(x)).collect(Collectors.toList());
    }

    private static List<List<HlaAllele>> cartesianAlleleProduct(final List<HlaAllele> first, final List<HlaAllele> second)
    {
        List<List<HlaAllele>> results = Lists.newArrayList();

        for(HlaAllele i : first)
        {
            for(HlaAllele j : second)
            {
                if(i != j)
                {
                    List<HlaAllele> pairing = Lists.newArrayList(i, j);
                    Collections.sort(pairing);
                    if(results.stream().anyMatch(x -> x.get(0).matches(pairing.get(0)) && x.get(1).matches(pairing.get(1))))
                        continue;

                    results.add(pairing);
                }
            }
        }

        return results;
    }

    private static List<HlaAllele> filterWithConfirmedProteins(final List<HlaAllele> alleles, List<HlaAllele> confirmedGroups)
    {
        return filterWithConfirmedGroups(alleles, confirmedGroups);
    }

    private static List<HlaAllele> filterWithConfirmedGroups(final List<HlaAllele> alleles, final List<HlaAllele> confirmedGroups)
    {
        // List<HlaAllele> alleleGroups = alleles.stream().map(x -> x.asAlleleGroup()).collect(Collectors.toList());
        //return filterWithConfirmed(alleleGroups, confirmedGroups);

        Map<String,List<HlaAllele>> map = Maps.newHashMap();

        map.put("A", confirmedGroups.stream().filter(x -> x.Gene.equals("A")).collect(Collectors.toList()));
        map.put("B", confirmedGroups.stream().filter(x -> x.Gene.equals("B")).collect(Collectors.toList()));
        map.put("C", confirmedGroups.stream().filter(x -> x.Gene.equals("C")).collect(Collectors.toList()));

        // return this.filter { map[it.gene]!!.size < 2 || map[it.gene]!!.contains(transform(it)) }

        List<HlaAllele> results = Lists.newArrayList();
        for(HlaAllele allele : alleles)
        {
            List<HlaAllele> geneAlleleList = map.get(allele.Gene);

            if(geneAlleleList.size() < 2 || contains(geneAlleleList, allele.asAlleleGroup()))
                results.add(allele);
        }

        return results;
    }

    private static List<HlaAllele> alleles(final List<HlaAlleleCoverage> coverage)
    {
        return coverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }
}
