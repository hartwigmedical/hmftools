package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.COMPLEX_PERMS_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.groupCoverage;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.proteinCoverage;
import static com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage.coverageAlleles;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;

public class HlaComplexBuilder
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    public HlaComplexBuilder(final LilacConfig config, final ReferenceData refData)
    {
        mConfig = config;
        mRefData = refData;
    }

    public List<HlaComplex> buildComplexes(final List<FragmentAlleles> referenceFragmentAlleles, final List<HlaAllele> candidateAlleles)
    {
        LL_LOGGER.info("Identifying uniquely identifiable groups and proteins [total,unique,shared,wild]");

        HlaComplexCoverage groupCoverage = groupCoverage(referenceFragmentAlleles, candidateAlleles);

        List<HlaAlleleCoverage> confirmedGroups = confirmUnique(groupCoverage, Lists.newArrayList());

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

        List<HlaAllele> confirmedGroupAlleles = coverageAlleles(confirmedGroups);
        List<HlaAllele> candidatesAfterConfirmedGroups = filterWithConfirmedGroups(candidateAlleles, confirmedGroupAlleles);
        HlaComplexCoverage proteinCoverage = proteinCoverage(referenceFragmentAlleles, candidatesAfterConfirmedGroups);
        List<HlaAlleleCoverage> confirmedProtein = confirmUnique(proteinCoverage,  confirmedGroupAlleles);
        List<HlaAlleleCoverage> discardedProtein = proteinCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !confirmedProtein.contains(x)).collect(Collectors.toList());
        Collections.sort(discardedProtein, Collections.reverseOrder());

        if(!confirmedProtein.isEmpty())
        {
            LL_LOGGER.info("  confirmed {} unique proteins: {}", confirmedProtein.size(), HlaAlleleCoverage.toString(confirmedProtein));
        }
        else
        {
            LL_LOGGER.info("  confirmed 0 unique proteins");
        }

        if (!discardedProtein.isEmpty())
        {
            LL_LOGGER.info("  found {} insufficiently unique proteins: {}", discardedProtein.size(), HlaAlleleCoverage.toString(discardedProtein));
        }

        List<HlaAllele> confirmedProteinAlleles = coverageAlleles(confirmedProtein);

        List<HlaAllele> candidatesAfterConfirmedProteins = filterWithConfirmedProteins(candidatesAfterConfirmedGroups, confirmedProteinAlleles);

        List<HlaComplex> aOnlyComplexes = buildComplexesByGene(GENE_A, confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);
        List<HlaComplex> bOnlyComplexes = buildComplexesByGene(GENE_B, confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);
        List<HlaComplex> cOnlyComplexes = buildComplexesByGene(GENE_C, confirmedGroupAlleles, confirmedProteinAlleles, candidatesAfterConfirmedProteins);

        List<HlaComplex> complexes;
        long simpleComplexCount = (long)aOnlyComplexes.size() * bOnlyComplexes.size() * cOnlyComplexes.size();

        if (simpleComplexCount > COMPLEX_PERMS_THRESHOLD || simpleComplexCount < 0)
        {
            LL_LOGGER.info("candidate permutations exceeds maximum complexity, complexes(A={} B={} C={})",
                    aOnlyComplexes.size(), bOnlyComplexes.size(), cOnlyComplexes.size());

            List<HlaAllele> aTopCandidates = rankedGroupCoverage(10, referenceFragmentAlleles, aOnlyComplexes);
            List<HlaAllele> bTopCandidates = rankedGroupCoverage(10, referenceFragmentAlleles, bOnlyComplexes);
            List<HlaAllele> cTopCandidates = rankedGroupCoverage(10, referenceFragmentAlleles, cOnlyComplexes);
            List<HlaAllele> topCandidates = Lists.newArrayList();
            topCandidates.addAll(aTopCandidates);
            topCandidates.addAll(bTopCandidates);
            topCandidates.addAll(cTopCandidates);

            List<HlaAllele> rejected = candidatesAfterConfirmedProteins.stream()
                    .filter(x -> !topCandidates.contains(x)).collect(Collectors.toList());

            LL_LOGGER.info("  discarding {} unlikely candidates: {}", rejected.size(), HlaAllele.toString(rejected));

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
        List<HlaComplex> a = buildComplexesByGene(GENE_A, confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> b = buildComplexesByGene(GENE_B, confirmedGroups, confirmedProteins, candidates);
        List<HlaComplex> c = buildComplexesByGene(GENE_C, confirmedGroups, confirmedProteins, candidates);
        return combineComplexes(combineComplexes(a, b), c);
    }

    public static List<HlaComplex> buildComplexesByGene(
            final String gene, final List<HlaAllele> unfilteredGroups,
            final List<HlaAllele> unfilteredProteins, final List<HlaAllele> unfilteredCandidates)
    {
        List<HlaAllele> confirmedGroups = takeN(unfilteredGroups.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList()), 2);
        List<HlaAllele> confirmedProteins = takeN(unfilteredProteins.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList()), 2);
        List<HlaAllele> candidates = unfilteredCandidates.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList());

        if (confirmedProteins.size() == 2)
            return Lists.newArrayList(new HlaComplex(confirmedProteins));

        if (confirmedProteins.size() == 1)
        {
            List<HlaAllele> confirmedProteinGroups = confirmedProteins.stream().map(x -> x.asAlleleGroup()).collect(Collectors.toList());
            List<HlaAllele> remainingGroups = confirmedGroups.stream().filter(x -> !confirmedProteinGroups.contains(x)).collect(Collectors.toList());

            List<HlaAllele> first = confirmedProteins;

            List<HlaAllele> second = remainingGroups.isEmpty() ?
                    candidates.stream().filter(x -> x != confirmedProteins.get(0)).collect(Collectors.toList()) :
                    candidates.stream().filter(x -> remainingGroups.contains(x.asAlleleGroup())).collect(Collectors.toList());

            List<HlaComplex> complexes = combineAlleles(first, second);
            if(!remainingGroups.isEmpty())
                return complexes;

            complexes.add(new HlaComplex(first));
            return complexes;

            // return if (remainingGroups.isEmpty()) combineAlleles(first, second) + HlaComplex(first) else combineAlleles(first, second)
        }

        if (confirmedGroups.size() == 2)
        {
            List<HlaAllele> first = candidates.stream().filter(x -> x.asAlleleGroup() == confirmedGroups.get(0)).collect(Collectors.toList());
            List<HlaAllele> second = candidates.stream().filter(x -> x.asAlleleGroup() == confirmedGroups.get(1)).collect(Collectors.toList());
            return combineAlleles(first, second);
        }

        if (confirmedGroups.size() == 1)
        {
            List<HlaAllele> first = candidates.stream().filter(x -> x.asAlleleGroup() == confirmedGroups.get(0)).collect(Collectors.toList());
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
        // first produce each unique combo pairing, then combine into a single complex
        List<List<HlaComplex>> intermediatePairs = cartesianComplexProduct(first, second);

        List<HlaComplex> complexes = Lists.newArrayList();

        for(List<HlaComplex> pairing : intermediatePairs)
        {
            List<HlaAllele> combinedAlleles = pairing.get(0).getAlleles().stream().collect(Collectors.toList());
            combinedAlleles.addAll(pairing.get(1).getAlleles());
            complexes.add(new HlaComplex(combinedAlleles));
        }

        return complexes;
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

                    // this check can be avoided given than A/B, or B/Cs are always being combined
                    //if(results.stream().anyMatch(x -> x.get(0) == pairing.get(0) && x.get(1) == pairing.get(1)))
                    //    continue;

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
        // make a list of all possible combinations of the alleles in each of the 2 lists
        List<List<HlaAllele>> results = Lists.newArrayList();

        for(HlaAllele i : first)
        {
            for(HlaAllele j : second)
            {
                if(i != j)
                {
                    List<HlaAllele> pairing = Lists.newArrayList(i, j);
                    Collections.sort(pairing);
                    if(results.stream().anyMatch(x -> x.get(0) == pairing.get(0) && x.get(1) == pairing.get(1)))
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
        Map<String,List<HlaAllele>> map = Maps.newHashMap();

        GENE_IDS.forEach(x -> map.put(x, confirmedGroups.stream().filter(y -> y.Gene.equals(x)).collect(Collectors.toList())));

        List<HlaAllele> results = Lists.newArrayList();
        for(HlaAllele allele : alleles)
        {
            List<HlaAllele> geneAlleleList = map.get(allele.Gene);

            if(geneAlleleList.size() < 2 || geneAlleleList.contains(allele.asAlleleGroup()))
                results.add(allele);
        }

        return results;
    }

    private List<HlaAlleleCoverage> confirmUnique(final HlaComplexCoverage complexCoverage, final List<HlaAllele> confirmedGroupAlleles)
    {
        List<HlaAlleleCoverage> unique = complexCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage >= mConfig.MinConfirmedUniqueCoverage).collect(Collectors.toList());

        Collections.sort(unique, Collections.reverseOrder());

        List<HlaAlleleCoverage> results = Lists.newArrayList();

        // take at most 2 alleles for each gene, and at most 1 unique protein if more than 1 unique group is provided
        for(String gene : GENE_IDS)
        {
            List<HlaAlleleCoverage> geneCoverage = unique.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

            if(geneCoverage.isEmpty())
                continue;

            final List<HlaAllele> geneGroupAlleles = confirmedGroupAlleles.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList());

            int geneCount = 0;
            for(HlaAlleleCoverage coverage : geneCoverage)
            {
                if(geneGroupAlleles.size() > 1)
                {
                    // how many added already from this protein's group
                    int matchedGroupCount = (int)results.stream()
                            .filter(x -> geneGroupAlleles.contains(x.Allele.asAlleleGroup())
                                    && x.Allele.asAlleleGroup().equals(coverage.Allele.asAlleleGroup())).count();

                    if(matchedGroupCount >= 1)
                        continue;
                }

                results.add(coverage);

                ++geneCount;

                if(geneCount >= 2)
                    break;
            }
        }

        Collections.sort(results, Collections.reverseOrder());
        return results;
    }

    private List<HlaAllele> rankedGroupCoverage(
            int take, final List<FragmentAlleles> fragmentAlleles, final List<HlaComplex> complexes)
    {
        List<HlaComplexCoverage> complexCoverages = complexes.stream()
                .map(x -> proteinCoverage(fragmentAlleles, x.getAlleles())).collect(Collectors.toList());

        HlaComplexCoverageRanking complexRanker = new HlaComplexCoverageRanking(0, mRefData);
        complexCoverages = complexRanker.rankCandidates(complexCoverages);

        // Collections.sort(complexCoverages, new HlaComplexCoverage.TotalCoverageSorter());
        List<HlaAllele> topRanked = Lists.newArrayList();

        for(HlaComplexCoverage coverage : complexCoverages)
        {
            coverage.getAlleleCoverage().stream().filter(x -> !topRanked.contains(x.Allele)).forEach(x -> topRanked.add(x.Allele));

            if(topRanked.size() >= take)
                break;
        }

        return topRanked;

        /* no longer reincluding common alleles if they ranked too low

        List<HlaAllele> topTakers = topRanked;

        List<HlaAllele> topRankedKeepers = topRanked.stream().filter(x -> mRefData.CommonAlleles.contains(x)).collect(Collectors.toList());

        List<HlaAllele> uniqueRanked = topTakers.stream().collect(Collectors.toList());
        topRankedKeepers.stream().filter(x -> !topTakers.contains(x)).forEach(x -> uniqueRanked.add(x));

        return uniqueRanked;
        */
    }

    private static List<HlaAllele> takeN(final List<HlaAllele> list, int n)
    {
        List<HlaAllele> newList = Lists.newArrayList();

        for(int i = 0; i < min(list.size(), n); ++i)
        {
            newList.add(list.get(i));
        }

        return newList;
    }

    /*
    private static List<HlaAlleleCoverage> takeN(final List<HlaAlleleCoverage> list, int n)
    {
        List<HlaAlleleCoverage> newList = Lists.newArrayList();

        for(int i = 0; i < min(list.size(), n); ++i)
        {
            newList.add(list.get(i));
        }

        return newList;
    }
    */

}
