package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.COMPLEX_PERMS_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_HI_CONF_UNIQUE_GROUP_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_LOW_CONF_GROUP_TOTAL_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.RARE_ALLELES_FREQ_CUTOFF;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.coverage.AlleleCoverage.coverageAlleles;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.filterUnsupportedWildcardFragments;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findUnsupportedWildcards;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findWildcardAlleles;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;

import org.apache.commons.lang3.tuple.Pair;

public class ComplexBuilder
{
    private final ReferenceData mRefData;

    private final List<HlaAllele> mUniqueGroupAlleles;
    private final List<HlaAllele> mConfirmedProteinAlleles;
    private final List<HlaAllele> mUniqueProteinAlleles;
    private final List<HlaAllele> mConfirmedRecoveredAlleles;
    private final int mGeneCount;

    private static final int DIPLOID_ALLELE_COUNT = 2;

    public ComplexBuilder(final ReferenceData refData, int geneCount)
    {
        mRefData = refData;

        mUniqueGroupAlleles = Lists.newArrayList();
        mUniqueProteinAlleles = Lists.newArrayList();
        mConfirmedProteinAlleles = Lists.newArrayList();
        mConfirmedRecoveredAlleles = Lists.newArrayList();
        mGeneCount = geneCount;
    }

    public List<HlaAllele> getUniqueProteinAlleles() { return mUniqueProteinAlleles; }
    public List<HlaAllele> getConfirmedRecoveredAlleles() { return mConfirmedRecoveredAlleles; }

    public void filterCandidates(
            final List<FragmentAlleles> refFragAlleles, final List<HlaAllele> candidateAlleles, final List<HlaAllele> recoveredAlleles)
    {
        LL_LOGGER.info("filtering candidates from fragAlleles({}) candidates({}) recovered({})",
                refFragAlleles.size(), candidateAlleles.size(), recoveredAlleles.size());

        ComplexCoverage groupCoverage = calcGroupCoverage(refFragAlleles, candidateAlleles);

        int totalFragCount = refFragAlleles.size();

        List<AlleleCoverage> uniqueGroups = findUnique(groupCoverage, Lists.newArrayList(), totalFragCount);

        double lowConfGroupMinUniqueFrags = totalFragCount * MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE;
        double lowConfGroupMinTotalFrags = totalFragCount * MIN_LOW_CONF_GROUP_TOTAL_COVERAGE;

        List<AlleleCoverage> lowConfGroups = Lists.newArrayList();
        List<AlleleCoverage> discardedGroups = Lists.newArrayList();

        for(AlleleCoverage alleleCoverage : groupCoverage.getAlleleCoverage())
        {
            if(uniqueGroups.contains(alleleCoverage))
                continue;

            if(alleleCoverage.UniqueCoverage >= lowConfGroupMinUniqueFrags || alleleCoverage.TotalCoverage >= lowConfGroupMinTotalFrags)
            {
                lowConfGroups.add(alleleCoverage);
            }
            else
            {
                discardedGroups.add(alleleCoverage);
            }
        }

        LL_LOGGER.info("  confirmed {} unique allele groups{}{}",
                uniqueGroups.size(), uniqueGroups.isEmpty() ? "" : ": ", AlleleCoverage.toString(uniqueGroups));

        LL_LOGGER.info("  found {} insufficiently unique allele groups{}{}",
                lowConfGroups.size(), lowConfGroups.isEmpty() ? "" : ": ", AlleleCoverage.toString(lowConfGroups));

        LL_LOGGER.debug("  discarded {} allele groups{}{}",
                discardedGroups.size(), discardedGroups.isEmpty() ? "" : ": ", AlleleCoverage.toString(discardedGroups));

        List<HlaAllele> uniqueGroupAlleles = coverageAlleles(uniqueGroups);

        final List<HlaAllele> candidatesAfterUniqueGroups = filterWithUniqueGroups(
                candidateAlleles, uniqueGroupAlleles, Sets.newHashSet(recoveredAlleles));

        // keep common alleles in insufficiently unique groups
        List<HlaAllele> topLowConfGroups = getTopLowConfGroups(uniqueGroups, lowConfGroups);
        List<HlaAllele> commonAllelesInLowConfGroups = mRefData.CommonAlleles.stream()
                .filter(x -> !candidatesAfterUniqueGroups.contains(x))
                .filter(x -> topLowConfGroups.contains(x.asAlleleGroup()))
                .collect(Collectors.toList());

        candidatesAfterUniqueGroups.addAll(commonAllelesInLowConfGroups);

        LL_LOGGER.debug("  keeping {} common allele(s) from insufficiently unique allele groups{}{}",
                commonAllelesInLowConfGroups.size(), commonAllelesInLowConfGroups.isEmpty() ? "" : ": ",
                HlaAllele.toString(commonAllelesInLowConfGroups));

        // keep common alleles from the same 2-digit group as rare alleles or alleles with wildcards
        List<HlaAllele> rareAlleles = candidatesAfterUniqueGroups.stream()
                .filter(x -> mRefData.getAlleleFrequencies().getAlleleFrequency(x) <= RARE_ALLELES_FREQ_CUTOFF)
                .toList();

        Set<HlaAllele> alleleGroupsOfDubiousAlleles = candidatesAfterUniqueGroups.stream()
                .filter(x -> x.hasWildcards() || rareAlleles.contains(x))
                .map(HlaAllele::asAlleleGroup)
                .collect(Collectors.toSet());

        List<HlaAllele> commonAllelesFromSameGroupAsDubiousAlleles = mRefData.CommonAlleles.stream()
                .filter(x -> !candidatesAfterUniqueGroups.contains(x))
                .filter(x -> alleleGroupsOfDubiousAlleles.contains(x.asAlleleGroup()))
                .collect(Collectors.toList());

        candidatesAfterUniqueGroups.addAll(commonAllelesFromSameGroupAsDubiousAlleles);

        LL_LOGGER.debug("  keeping {} common allele(s) in the same 2-digit group as wildcard or rare allele candidates{}{}",
                commonAllelesFromSameGroupAsDubiousAlleles.size(),
                commonAllelesFromSameGroupAsDubiousAlleles.isEmpty() ? "" : ": ",
                HlaAllele.toString(commonAllelesFromSameGroupAsDubiousAlleles));

        // keep known stop-loss alleles
        List<HlaAllele> allelesWithStopLossIndel = recoveredAlleles.stream()
                .filter(mRefData.KnownStopLossIndelAlleles::containsValue)
                .filter(x -> !candidatesAfterUniqueGroups.contains(x))
                .collect(Collectors.toList());

        candidatesAfterUniqueGroups.addAll(allelesWithStopLossIndel);

        LL_LOGGER.debug("  keeping {} allele(s) with stop loss indel{}{}",
                allelesWithStopLossIndel.size(), allelesWithStopLossIndel.isEmpty() ? "" : ": ",
                HlaAllele.toString(allelesWithStopLossIndel));

        recoveredAlleles.stream().filter(candidatesAfterUniqueGroups::contains).forEach(mConfirmedRecoveredAlleles::add);

        if(!mConfirmedRecoveredAlleles.isEmpty())
        {
            Collections.sort(mConfirmedRecoveredAlleles);
            LL_LOGGER.info("  keeping {} recovered alleles from insufficiently unique groups: {}",
                    mConfirmedRecoveredAlleles.size(), HlaAllele.toString(mConfirmedRecoveredAlleles));
        }
        else if(!recoveredAlleles.isEmpty())
        {
            LL_LOGGER.info("  no recovered alleles kept from insufficiently unique groups");
        }

        ComplexCoverage proteinCoverage = calcProteinCoverage(refFragAlleles, candidatesAfterUniqueGroups);

        // find uniquely supported protein alleles but don't allow recovered alleles to be in the unique protein set
        List<AlleleCoverage> uniqueProteins = findUnique(proteinCoverage, uniqueGroupAlleles, totalFragCount).stream()
                .filter(x -> !recoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        List<AlleleCoverage> insufficientlyUniqueProteins = proteinCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage > 0 && !uniqueProteins.contains(x))
                .sorted(Collections.reverseOrder())
                .collect(Collectors.toList());

        LL_LOGGER.info("  confirmed {} unique proteins{}{}",
                uniqueProteins.size(), uniqueProteins.isEmpty() ? "" : ": ", AlleleCoverage.toString(uniqueProteins));

        LL_LOGGER.info("  found {} insufficiently unique proteins{}{}",
                insufficientlyUniqueProteins.size(), insufficientlyUniqueProteins.isEmpty() ? "" : ": ",
                AlleleCoverage.toString(insufficientlyUniqueProteins));

        // unique protein filtering is no longer applied

        mUniqueProteinAlleles.addAll(filterWithUniqueProteins(candidatesAfterUniqueGroups, mConfirmedProteinAlleles));

        mUniqueGroupAlleles.addAll(uniqueGroupAlleles);
    }

    public List<HlaComplex> buildComplexes(final List<FragmentAlleles> refFragAlleles, final Set<HlaAllele> recoveredAlleles)
    {
        // filter out any wildcards
        Set<HlaAllele> wildcardAlleles = findWildcardAlleles(refFragAlleles);
        List<HlaAllele> unsupportedWildcards = findUnsupportedWildcards(refFragAlleles, wildcardAlleles);
        filterUnsupportedWildcardFragments(refFragAlleles, unsupportedWildcards);

        unsupportedWildcards.forEach(mConfirmedProteinAlleles::remove);
        unsupportedWildcards.forEach(mUniqueProteinAlleles::remove);

        Map<HlaGene, List<HlaComplex>> geneOnlyComplexes = Maps.newHashMap();
        for(HlaGene gene : GENE_CACHE.GeneNames)
        {
            if(gene.isPseudo())
                continue;

            geneOnlyComplexes.put(gene, buildComplexesByGene(gene, mUniqueGroupAlleles, mUniqueProteinAlleles));
        }

        List<HlaComplex> complexes;
        long simpleComplexCount = geneOnlyComplexes.values().stream().mapToLong(List::size).reduce(1L, (acc, x) -> acc * x);

        if(simpleComplexCount > COMPLEX_PERMS_THRESHOLD || simpleComplexCount < 0)
        {
            // common alleles which satisfy the filtering by unique groups will be kept regardless of any ranking
            List<HlaAllele> commonAlleles = mRefData.CommonAlleles.stream()
                    .filter(mUniqueProteinAlleles::contains)
                    .toList();

            Map<HlaGene, Integer> geneOnlyComplexeSizes =
                    geneOnlyComplexes.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, x -> x.getValue().size()));
            LL_LOGGER.info("candidate permutations exceeds threshold, candidates({}) common({})",
                    geneOnlyComplexeSizes, commonAlleles.size());

            Map<HlaGene, List<HlaAllele>> geneTopCandidates = Maps.newHashMap();
            for(Map.Entry<HlaGene, List<HlaComplex>> entry : geneOnlyComplexes.entrySet())
                geneTopCandidates.put(entry.getKey(), rankedGroupCoverage(10, refFragAlleles, entry.getValue(), recoveredAlleles));

            List<HlaAllele> topCandidates = geneTopCandidates.values().stream()
                    .flatMap(Collection::stream)
                    .collect(Collectors.toCollection(Lists::newArrayList));

            // ensure any common alleles, unfiltered so far, are kept regardless of ranking
            commonAlleles.stream().filter(x -> !topCandidates.contains(x)).forEach(topCandidates::add);

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

    private static ComplexCoverage calcGroupCoverage(final List<FragmentAlleles> fragAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(fragAlleles, alleles);
        return ComplexCoverage.create(AlleleCoverage.groupCoverage(filteredFragments));
    }

    public static ComplexCoverage calcProteinCoverage(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        List<FragmentAlleles> filteredFragments = FragmentAlleles.filter(fragmentAlleles, alleles);
        return ComplexCoverage.create(AlleleCoverage.proteinCoverage(filteredFragments));
    }

    private static List<HlaComplex> buildAlleleComplexes(
            final Collection<HlaAllele> confirmedGroups, final Collection<HlaAllele> candidates)
    {
        return GENE_CACHE.GeneNames.stream()
                .filter(gene -> !gene.isPseudo())
                .map(gene -> buildComplexesByGene(gene, confirmedGroups, candidates))
                .reduce(ComplexBuilder::combineComplexes)
                .orElse(null);
    }

    public static List<HlaComplex> buildComplexesByGene(
            final HlaGene gene, final Collection<HlaAllele> unfilteredGroups, final Collection<HlaAllele> unfilteredCandidates)
    {
        List<HlaAllele> confirmedGroups = takeN(unfilteredGroups.stream().filter(x -> x.Gene == gene).collect(Collectors.toList()), 2);
        List<HlaAllele> candidates = unfilteredCandidates.stream().filter(x -> x.Gene == gene).collect(Collectors.toList());
        Set<HlaAllele> candidatesSet = Sets.newHashSet(candidates);

        if(confirmedGroups.size() == 2)
        {
            Set<HlaAllele> first = candidates.stream()
                    .filter(x -> x.asAlleleGroup() == confirmedGroups.get(0))
                    .collect(Collectors.toCollection(Sets::newHashSet));
            Set<HlaAllele> second = candidates.stream()
                    .filter(x -> x.asAlleleGroup() == confirmedGroups.get(1))
                    .collect(Collectors.toCollection(Sets::newHashSet));
            return combineAlleles(first, second);
        }

        if(confirmedGroups.size() == 1)
        {
            Set<HlaAllele> first = candidates.stream()
                    .filter(x -> x.asAlleleGroup() == confirmedGroups.get(0))
                    .collect(Collectors.toCollection(Sets::newHashSet));
            Set<HlaAllele> second = candidatesSet;

            List<HlaComplex> complexes = first.stream().map(x -> new HlaComplex(Lists.newArrayList(x))).collect(Collectors.toList());
            complexes.addAll(combineAlleles(first, second));
            return complexes;
        }

        List<HlaComplex> complexes = candidates.stream().map(x -> new HlaComplex(Lists.newArrayList(x))).collect(Collectors.toList());
        complexes.addAll(combineAlleles(candidatesSet, candidatesSet));
        return complexes;
    }

    private static List<HlaComplex> combineComplexes(final Iterable<HlaComplex> first, final Iterable<HlaComplex> second)
    {
        // first produce each unique combo pairing, then combine into a single complex
        List<List<HlaComplex>> intermediatePairs = cartesianComplexProduct(first, second);

        List<HlaComplex> complexes = Lists.newArrayList();

        for(List<HlaComplex> pairing : intermediatePairs)
        {
            List<HlaAllele> combinedAlleles = Lists.newArrayList(pairing.get(0).Alleles);
            combinedAlleles.addAll(pairing.get(1).Alleles);
            complexes.add(new HlaComplex(combinedAlleles));
        }

        return complexes;
    }

    private static List<List<HlaComplex>> cartesianComplexProduct(final Iterable<HlaComplex> first, final Iterable<HlaComplex> second)
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

    private static List<HlaComplex> combineAlleles(final Set<HlaAllele> first, final Set<HlaAllele> second)
    {
        Set<Pair<HlaAllele, HlaAllele>> allelePairs = cartesianAlleleProduct(first, second);
        return allelePairs.stream().map(x -> new HlaComplex(Lists.newArrayList(x.getLeft(), x.getRight()))).collect(Collectors.toList());
    }

    private static Set<Pair<HlaAllele, HlaAllele>> cartesianAlleleProduct(final Set<HlaAllele> first, final Set<HlaAllele> second)
    {
        // make a list of all possible combinations of the alleles in each of the 2 lists
        Set<Pair<HlaAllele, HlaAllele>> results = Sets.newHashSet();

        for(HlaAllele allele1 : first)
        {
            for(HlaAllele allele2 : second)
            {
                if(allele1 != allele2)
                {
                    Pair<HlaAllele, HlaAllele> pairing;
                    if(allele1.compareTo(allele2) <= 0)
                        pairing = Pair.of(allele1, allele2);
                    else
                        pairing = Pair.of(allele2, allele1);

                    results.add(pairing);
                }
            }
        }

        return results;
    }

    private static List<HlaAllele> filterWithUniqueProteins(final List<HlaAllele> alleles, final Collection<HlaAllele> confirmedGroups)
    {
        return filterWithUniqueGroups(alleles, confirmedGroups, Sets.newHashSet());
    }

    private static List<HlaAllele> filterWithUniqueGroups(
            final Iterable<HlaAllele> alleles, final Collection<HlaAllele> confirmedGroups, final Set<HlaAllele> recoveredAlleles)
    {
        Map<HlaGene, List<HlaAllele>> map = Maps.newHashMap();

        GENE_CACHE.GeneNames.forEach(x -> map.put(x, confirmedGroups.stream().filter(y -> y.Gene == x).collect(Collectors.toList())));

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

    private static List<HlaAllele> getTopLowConfGroups(
            final Iterable<AlleleCoverage> uniqueGroups, final Iterable<AlleleCoverage> lowConfGroups)
    {
        Map<HlaGene, Integer> groupCountsPerGene = Maps.newHashMap();
        for(AlleleCoverage alleleCoverage : uniqueGroups)
            groupCountsPerGene.merge(alleleCoverage.Allele.Gene, 1, Integer::sum);

        List<AlleleCoverage> topLowConfGroups = Lists.newArrayList();
        for(AlleleCoverage alleleCoverage : lowConfGroups)
        {
            HlaGene gene = alleleCoverage.Allele.Gene;
            groupCountsPerGene.putIfAbsent(gene, 0);

            if(groupCountsPerGene.get(gene) < DIPLOID_ALLELE_COUNT)
                topLowConfGroups.add(alleleCoverage);

            groupCountsPerGene.merge(gene, 1, Integer::sum);
        }

        return coverageAlleles(topLowConfGroups);
    }

    private static double requiredUniqueGroupCoverage(double totalCoverage, boolean isGroup)
    {
        return isGroup ? totalCoverage * MIN_HI_CONF_UNIQUE_GROUP_COVERAGE : totalCoverage * MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE;
    }

    private static List<AlleleCoverage> findUnique(
            final ComplexCoverage complexCoverage, final Collection<HlaAllele> confirmedGroupAlleles, int totalFragCount)
    {
        List<AlleleCoverage> unique = complexCoverage.getAlleleCoverage().stream()
                .filter(x -> x.UniqueCoverage >= requiredUniqueGroupCoverage(totalFragCount, confirmedGroupAlleles.isEmpty()))
                .sorted(Collections.reverseOrder())
                .toList();

        List<AlleleCoverage> results = Lists.newArrayList();

        // take at most 2 alleles for each gene, and at most 1 unique protein if more than 1 unique group is provided
        for(HlaGene gene : GENE_CACHE.GeneNames)
        {
            List<AlleleCoverage> geneCoverage = unique.stream().filter(x -> x.Allele.Gene == gene).toList();

            if(geneCoverage.isEmpty())
                continue;

            final List<HlaAllele> geneGroupAlleles = confirmedGroupAlleles.stream()
                    .filter(x -> x.Gene == gene)
                    .toList();

            int geneCount = 0;
            for(AlleleCoverage coverage : geneCoverage)
            {
                if(geneGroupAlleles.size() > 1)
                {
                    // how many added already from this protein's group
                    int matchedGroupCount = (int) results.stream()
                            .filter(x -> geneGroupAlleles.contains(x.Allele.asAlleleGroup())
                                    && x.Allele.asAlleleGroup().equals(coverage.Allele.asAlleleGroup())).count();

                    if(matchedGroupCount >= 1)
                        continue;
                }

                results.add(coverage);

                ++geneCount;

                if(geneCount >= DIPLOID_ALLELE_COUNT)
                    break;
            }
        }

        Collections.sort(results, Collections.reverseOrder());
        return results;
    }

    private List<HlaAllele> rankedGroupCoverage(
            int take, final List<FragmentAlleles> fragAlleles, final Collection<HlaComplex> complexes,
            final Set<HlaAllele> recoveredAlleles)
    {
        List<ComplexCoverage> complexCoverages = complexes.stream()
                .map(x -> calcProteinCoverage(fragAlleles, x.Alleles)).collect(Collectors.toList());

        ComplexCoverageRanking complexRanker = new ComplexCoverageRanking(0, mRefData, mGeneCount);
        complexCoverages = complexRanker.rankCandidates(complexCoverages, recoveredAlleles, Lists.newArrayList());

        // take the top N alleles but no more than 5 that pair with something in the top 10
        Map<HlaAllele, Integer> pairingCount = Maps.newHashMap();

        List<HlaAllele> topRanked = Lists.newArrayList();

        for(ComplexCoverage coverage : complexCoverages)
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
