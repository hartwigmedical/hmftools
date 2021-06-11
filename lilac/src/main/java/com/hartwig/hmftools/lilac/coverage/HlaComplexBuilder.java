package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.COMPLEX_PERMS_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_CONF_UNIQUE_GROUP_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_CONF_UNIQUE_PROTEIN_COVERAGE;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.filterUnsupportedWildcardFragments;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findUnsupportedWildcards;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findWildcardAlleles;
import static com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage.coverageAlleles;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaComplexBuilder
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    private final List<HlaAllele> mUniqueGroupAlleles;
    private final List<HlaAllele> mConfirmedProteinAlleles;
    private final List<HlaAllele> mUniqueProteinAlleles;
    private final List<HlaAllele> mConfirmedRecoveredAlleles;

    public HlaComplexBuilder(final LilacConfig config, final ReferenceData refData)
    {
        mConfig = config;
        mRefData = refData;

        mUniqueGroupAlleles = Lists.newArrayList();
        mUniqueProteinAlleles = Lists.newArrayList();
        mConfirmedProteinAlleles = Lists.newArrayList();
        mConfirmedRecoveredAlleles = Lists.newArrayList();
    }

    public List<HlaAllele> getUniqueProteinAlleles() { return mUniqueProteinAlleles; }
    public List<HlaAllele> getConfirmedRecoveredAlleles() { return mConfirmedRecoveredAlleles; }

    public void filterCandidates(
            final List<FragmentAlleles> refFragAlleles, final List<HlaAllele> candidateAlleles, final List<HlaAllele> recoveredAlleles)
    {
        LL_LOGGER.info("filtering candidates from fragAlleles({}) candidates({}) recovered({})",
                refFragAlleles.size(), candidateAlleles.size(), recoveredAlleles.size());

        HlaComplexCoverage groupCoverage = calcGroupCoverage(refFragAlleles, candidateAlleles);

        int totalFragCount = refFragAlleles.size();
        List<HlaAlleleCoverage> uniqueGroups = findUnique(groupCoverage, Lists.newArrayList(), totalFragCount);

        List<HlaAlleleCoverage> discardedGroups = groupCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !uniqueGroups.contains(x)).collect(Collectors.toList());

        Collections.sort(discardedGroups, Collections.reverseOrder());

        if(!uniqueGroups.isEmpty())
        {
            LL_LOGGER.info("  confirmed {} unique groups: {}", uniqueGroups.size(), HlaAlleleCoverage.toString(uniqueGroups));
        }
        else
        {
            LL_LOGGER.info("  confirmed 0 unique groups");
        }

        if(!discardedGroups.isEmpty())
        {
            LL_LOGGER.info("  found {} insufficiently unique groups: {}",
                    discardedGroups.size(), HlaAlleleCoverage.toString(discardedGroups));
        }

        List<HlaAllele> uniqueGroupAlleles = coverageAlleles(uniqueGroups);

        final List<HlaAllele> candidatesAfterUniqueGroups = filterWithUniqueGroups(candidateAlleles, uniqueGroupAlleles, recoveredAlleles);

        // ensure common alleles in unique groups are kept
        List<HlaAllele> discardedUniqueGroupAlleles = coverageAlleles(discardedGroups);
        List<HlaAllele> commonAllelesInDiscardedUniqueGroups = mRefData.CommonAlleles.stream()
                .filter(x -> !candidatesAfterUniqueGroups.contains(x))
                .filter(x -> discardedUniqueGroupAlleles.contains(x.asAlleleGroup()))
                .collect(Collectors.toList());

        candidatesAfterUniqueGroups.addAll(commonAllelesInDiscardedUniqueGroups);

        // ensure known stop-loss alleles are kept
        recoveredAlleles.stream()
                .filter(x -> mRefData.KnownStopLossIndelAlleles.containsValue(x))
                .filter(x -> !candidatesAfterUniqueGroups.contains(x))
                .forEach(x -> candidatesAfterUniqueGroups.add(x));

        mConfirmedRecoveredAlleles.addAll(recoveredAlleles.stream()
                .filter(x -> candidatesAfterUniqueGroups.contains(x)).collect(Collectors.toList()));

        if(!mConfirmedRecoveredAlleles.isEmpty())
        {
            Collections.sort(mConfirmedRecoveredAlleles);
            LL_LOGGER.info("  keeping {} recovered alleles from unique groups: {}",
                    mConfirmedRecoveredAlleles.size(), HlaAllele.toString(mConfirmedRecoveredAlleles));
        }
        else if(!recoveredAlleles.isEmpty())
        {
            LL_LOGGER.info("  no recovered alleles kept from unique groups");
        }

        HlaComplexCoverage proteinCoverage = calcProteinCoverage(refFragAlleles, candidatesAfterUniqueGroups);

        // find uniquely supported protein alleles but don't allow recovered alleles to be in the unique protein set
        List<HlaAlleleCoverage> uniqueProteins = findUnique(proteinCoverage, uniqueGroupAlleles, totalFragCount).stream()
                .filter(x -> !recoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        List<HlaAlleleCoverage> discardedProtein = proteinCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !uniqueProteins.contains(x)).collect(Collectors.toList());
        Collections.sort(discardedProtein, Collections.reverseOrder());

        if(!uniqueProteins.isEmpty())
        {
            LL_LOGGER.info("  confirmed {} unique proteins: {}", uniqueProteins.size(), HlaAlleleCoverage.toString(uniqueProteins));
        }
        else
        {
            LL_LOGGER.info("  confirmed 0 unique proteins");
        }

        if(!discardedProtein.isEmpty())
        {
            LL_LOGGER.info("  found {} insufficiently unique proteins: {}", discardedProtein.size(), HlaAlleleCoverage.toString(discardedProtein));
        }

        // unique protein filtering is no longer applied
        // mConfirmedProteinAlleles.addAall(coverageAlleles(uniqueProteins));

        mUniqueProteinAlleles.addAll(filterWithUniqueProteins(candidatesAfterUniqueGroups, mConfirmedProteinAlleles));

        mUniqueGroupAlleles.addAll(uniqueGroupAlleles);
    }

    public List<HlaComplex> buildComplexes(final List<FragmentAlleles> refFragAlleles, final List<HlaAllele> recoveredAlleles)
    {
        // filter out any wildcards
        Set<HlaAllele> wildcardAlleles = findWildcardAlleles(refFragAlleles);
        List<HlaAllele> unsupportedWildcards = findUnsupportedWildcards(refFragAlleles, wildcardAlleles);
        filterUnsupportedWildcardFragments(refFragAlleles, unsupportedWildcards);

        unsupportedWildcards.forEach(x -> mConfirmedProteinAlleles.remove(x));
        unsupportedWildcards.forEach(x -> mUniqueProteinAlleles.remove(x));

        List<HlaComplex> aOnlyComplexes = buildComplexesByGene(GENE_A, mUniqueGroupAlleles, mUniqueProteinAlleles);
        List<HlaComplex> bOnlyComplexes = buildComplexesByGene(GENE_B, mUniqueGroupAlleles, mUniqueProteinAlleles);
        List<HlaComplex> cOnlyComplexes = buildComplexesByGene(GENE_C, mUniqueGroupAlleles, mUniqueProteinAlleles);

        List<HlaComplex> complexes;
        long simpleComplexCount = (long)aOnlyComplexes.size() * bOnlyComplexes.size() * cOnlyComplexes.size();

        if (simpleComplexCount > COMPLEX_PERMS_THRESHOLD || simpleComplexCount < 0)
        {
            // common alleles which satisfy the filtering by unique groups will be kept regardless of any ranking
            List<HlaAllele> commonAlleles = mRefData.CommonAlleles.stream()
                    .filter(x -> mUniqueProteinAlleles.contains(x))
                    .collect(Collectors.toList());

            LL_LOGGER.info("candidate permutations exceeds threshold, candidates(A={} B={} C={}) common({})",
                    aOnlyComplexes.size(), bOnlyComplexes.size(), cOnlyComplexes.size(), commonAlleles.size());

            List<HlaAllele> aTopCandidates = rankedGroupCoverage(10, refFragAlleles, aOnlyComplexes, recoveredAlleles);
            List<HlaAllele> bTopCandidates = rankedGroupCoverage(10, refFragAlleles, bOnlyComplexes, recoveredAlleles);
            List<HlaAllele> cTopCandidates = rankedGroupCoverage(10, refFragAlleles, cOnlyComplexes, recoveredAlleles);
            List<HlaAllele> topCandidates = Lists.newArrayList();
            topCandidates.addAll(aTopCandidates);
            topCandidates.addAll(bTopCandidates);
            topCandidates.addAll(cTopCandidates);

            // ensure any common alleles, unfiltered so far, are kept regardless of ranking
            commonAlleles.stream().filter(x -> !topCandidates.contains(x)).forEach(x -> topCandidates.add(x));

            List<HlaAllele> rejected = mUniqueProteinAlleles.stream()
                    .filter(x -> !topCandidates.contains(x)).collect(Collectors.toList());

            LL_LOGGER.info("  discarding {} unlikely candidates: {}", rejected.size(), HlaAllele.toString(rejected));

            complexes = buildAlleleComplexes(mUniqueGroupAlleles, topCandidates);
        }
        else
        {
            complexes = buildAlleleComplexes(mUniqueGroupAlleles, mUniqueProteinAlleles);
        }

        return complexes;
    }

    private static HlaComplexCoverage calcGroupCoverage(final List<FragmentAlleles> fragAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(fragAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(filteredFragments));
    }

    public static HlaComplexCoverage calcProteinCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(fragmentAlleles, alleles);
        return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments));
    }

    private static List<HlaComplex> buildAlleleComplexes(final List<HlaAllele> confirmedGroups, final List<HlaAllele> candidates)
    {
        List<HlaComplex> a = buildComplexesByGene(GENE_A, confirmedGroups, candidates);
        List<HlaComplex> b = buildComplexesByGene(GENE_B, confirmedGroups, candidates);
        List<HlaComplex> c = buildComplexesByGene(GENE_C, confirmedGroups, candidates);
        return combineComplexes(combineComplexes(a, b), c);
    }

    public static List<HlaComplex> buildComplexesByGene(
            final String gene, final List<HlaAllele> unfilteredGroups, final List<HlaAllele> unfilteredCandidates)
    {
        List<HlaAllele> confirmedGroups = takeN(unfilteredGroups.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList()), 2);
        List<HlaAllele> candidates = unfilteredCandidates.stream().filter(x -> x.Gene.equals(gene)).collect(Collectors.toList());

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
            List<HlaAllele> combinedAlleles = pairing.get(0).Alleles.stream().collect(Collectors.toList());
            combinedAlleles.addAll(pairing.get(1).Alleles);
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

    private static List<HlaAllele> filterWithUniqueProteins(final List<HlaAllele> alleles, List<HlaAllele> confirmedGroups)
    {
        return filterWithUniqueGroups(alleles, confirmedGroups, Lists.newArrayList());
    }

    private static List<HlaAllele> filterWithUniqueGroups(
            final List<HlaAllele> alleles, final List<HlaAllele> confirmedGroups, final List<HlaAllele> recoveredAlleles)
    {
        Map<String,List<HlaAllele>> map = Maps.newHashMap();

        GENE_IDS.forEach(x -> map.put(x, confirmedGroups.stream().filter(y -> y.Gene.equals(x)).collect(Collectors.toList())));

        List<HlaAllele> results = Lists.newArrayList();
        for(HlaAllele allele : alleles)
        {
            List<HlaAllele> geneAlleleList = map.get(allele.Gene);

            // discard recovered alleles which aren't in a unique group
            if(recoveredAlleles.contains(allele) && !geneAlleleList.contains(allele.asAlleleGroup()))
                continue;

            if(geneAlleleList.size() < 2 || geneAlleleList.contains(allele.asAlleleGroup()))
                results.add(allele);
        }

        return results;
    }

    private static double requiredUniqueGroupCoverage(double totalCoverage, boolean isGroup)
    {
        return isGroup ? totalCoverage * MIN_CONF_UNIQUE_GROUP_COVERAGE : totalCoverage * MIN_CONF_UNIQUE_PROTEIN_COVERAGE;
    }

    private List<HlaAlleleCoverage> findUnique(
            final HlaComplexCoverage complexCoverage, final List<HlaAllele> confirmedGroupAlleles, int totalFragCount)
    {
        List<HlaAlleleCoverage> unique = complexCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage >= requiredUniqueGroupCoverage(totalFragCount, confirmedGroupAlleles.isEmpty()))
                .collect(Collectors.toList());

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
            int take, final List<FragmentAlleles> fragAlleles, final List<HlaComplex> complexes, final List<HlaAllele> recoveredAlleles)
    {
        List<HlaComplexCoverage> complexCoverages = complexes.stream()
                .map(x -> calcProteinCoverage(fragAlleles, x.Alleles)).collect(Collectors.toList());

        HlaComplexCoverageRanking complexRanker = new HlaComplexCoverageRanking(0, mRefData);
        complexCoverages = complexRanker.rankCandidates(complexCoverages, recoveredAlleles, Lists.newArrayList());

        // take the top N alleles but no more than 5 that pair with something in the top 10
        Map<HlaAllele,Integer> pairingCount = Maps.newHashMap();

        List<HlaAllele> topRanked = Lists.newArrayList();

        for(HlaComplexCoverage coverage : complexCoverages)
        {
            HlaAllele allele1 = coverage.getAlleles().get(0);

            if(coverage.getAlleles().size() == 1)
            {
                if(!topRanked.contains(allele1))
                    topRanked.add(allele1);
            }
            else
            {
                HlaAllele allele2 = coverage.getAlleles().get(1);
                Integer count1 = pairingCount.get(allele1);
                Integer count2 = pairingCount.get(allele2);

                if(count1 != null && count1 >= 5 && count2 == null)
                    continue;

                if(count2 != null && count2 >= 5 && count1 == null)
                    continue;

                if(count1 == null)
                {
                    pairingCount.put(allele1, 1);

                    if(!topRanked.contains(allele1))
                        topRanked.add(allele1);
                }
                else
                {
                    pairingCount.put(allele1, count1 + 1);
                }

                if(count2 == null)
                {
                    pairingCount.put(allele2, 1);

                    if(!topRanked.contains(allele2))
                        topRanked.add(allele2);
                }
                else
                {
                    pairingCount.put(allele2, count2 + 1);
                }
            }

            if(topRanked.size() >= take)
                break;
        }

        return topRanked;
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

}
