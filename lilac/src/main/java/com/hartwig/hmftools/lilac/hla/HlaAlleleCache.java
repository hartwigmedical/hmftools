package com.hartwig.hmftools.lilac.hla;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.math3.util.Pair;

public class HlaAlleleCache
{
    // ensure the same HlaAlle is not created twice, allowing for object matching instead of string matching

    // map of group alleles to protein/4-digit alleles to list of synonymous/6-digit alleles
    private final Map<HlaAllele,Map<HlaAllele,List<HlaAllele>>> mAlleleMap;

    public HlaAlleleCache()
    {
        mAlleleMap = Maps.newHashMap();
    }

    public int groupCount() { return mAlleleMap.size(); }
    public int fourDigitCount() { return mAlleleMap.values().stream().mapToInt(x -> x.size()).sum(); }
    public int alleleCount()
    {
        return mAlleleMap.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
    }

    public HlaAllele request(final String alleleStr)
    {
        HlaAllele tmpAllele = HlaAllele.fromString(alleleStr);

        // first check if the allele has already been created
        Pair<HlaAllele,Map<HlaAllele,List<HlaAllele>>> group = getOrCreateGroup(tmpAllele.asAlleleGroup());
        Pair<HlaAllele,List<HlaAllele>> fourDigit = getOrCreateFourDigitList(group.getFirst(), group.getSecond(), tmpAllele.asFourDigit());
        HlaAllele match = fourDigit.getSecond().stream().filter(x -> x == tmpAllele).findFirst().orElse(null);

        if(match != null)
            return match;

        HlaAllele newAllele = new HlaAllele(
                tmpAllele.Gene, tmpAllele.AlleleGroup, tmpAllele.Protein, tmpAllele.Synonymous, tmpAllele.SynonymousNonCoding,
                fourDigit.getFirst(), group.getFirst());

        fourDigit.getSecond().add(newAllele);
        return newAllele;
    }

    public HlaAllele requestFourDigit(final String alleleStr)
    {
        HlaAllele tmpAllele = HlaAllele.fromString(alleleStr).asFourDigit();

        Pair<HlaAllele,Map<HlaAllele,List<HlaAllele>>> group = getOrCreateGroup(tmpAllele.asAlleleGroup());

        HlaAllele match = group.getSecond().keySet().stream().filter(x -> x.equals(tmpAllele)).findFirst().orElse(null);
        if(match != null)
            return match;

        HlaAllele newAllele = new HlaAllele(
                tmpAllele.Gene, tmpAllele.AlleleGroup, tmpAllele.Protein, tmpAllele.Synonymous, tmpAllele.SynonymousNonCoding,
                null, group.getKey());

        group.getSecond().put(newAllele, Lists.newArrayList());
        return newAllele;
    }

    public HlaAllele requestGroup(final HlaAllele requested)
    {
        // ensure in group format
        HlaAllele allele = requested.asAlleleGroup().matches(requested) ? requested : requested.asAlleleGroup();
        Pair<HlaAllele,Map<HlaAllele,List<HlaAllele>>> group = getOrCreateGroup(allele);
        return group.getFirst();
    }

    private Pair<HlaAllele,Map<HlaAllele,List<HlaAllele>>> getOrCreateGroup(final HlaAllele allele)
    {
        HlaAllele groupAllele = allele.asAlleleGroup();

        Map.Entry<HlaAllele,Map<HlaAllele,List<HlaAllele>>> entry = mAlleleMap.entrySet().stream()
                .filter(x -> x.getKey().equals(groupAllele)).findFirst().orElse(null);

        if(entry != null)
            return Pair.create(entry.getKey(), entry.getValue());

        Map<HlaAllele,List<HlaAllele>> groupMap = Maps.newHashMap();

        HlaAllele newGroup = new HlaAllele(
                allele.Gene, allele.AlleleGroup, "", "", "", null, null);

        mAlleleMap.put(newGroup, groupMap);
        return Pair.create(newGroup, groupMap);
    }

    private static Pair<HlaAllele,List<HlaAllele>> getOrCreateFourDigitList(
            final HlaAllele groupAllele, final Map<HlaAllele,List<HlaAllele>> groupMap, final HlaAllele allele)
    {
        HlaAllele fourDigitAllele = allele.asFourDigit();

        Map.Entry<HlaAllele,List<HlaAllele>> entry = groupMap.entrySet().stream()
                .filter(x -> x.getKey().equals(fourDigitAllele)).findFirst().orElse(null);

        if(entry != null)
            return Pair.create(entry.getKey(), entry.getValue());

        List<HlaAllele> alleleList = Lists.newArrayList();

        HlaAllele newFourDigit = new HlaAllele(
                allele.Gene, allele.AlleleGroup, allele.Protein, "", "", null, groupAllele);

        groupMap.put(newFourDigit, alleleList);
        return Pair.create(newFourDigit, alleleList);
    }

    public void rebuildProteinAlleles(final List<HlaAllele> list)
    {
        List<HlaAllele> cachedList = list.stream().map(x -> requestFourDigit(x.asFourDigit().toString())).collect(Collectors.toList());
        list.clear();
        list.addAll(cachedList);
    }

    @VisibleForTesting
    public List<HlaAllele> getAlleles()
    {
        List<HlaAllele> allAlleles = Lists.newArrayList();
        mAlleleMap.values().stream().forEach(x -> x.values().stream().forEach(y -> allAlleles.addAll(y)));
        return allAlleles;
    }

    public HlaAllele findFourDigitAllele(final String alleleStr)
    {
        // assumes four-digit
        for(Map<HlaAllele,List<HlaAllele>> proteinEntry : mAlleleMap.values())
        {
            HlaAllele matchedAllele = proteinEntry.keySet().stream().filter(x -> x.matches(alleleStr)).findFirst().orElse(null);
            if(matchedAllele != null)
                return matchedAllele;
        }

        return null;
    }

    public HlaAllele findAllele(final String alleleStr)
    {
        // assumes four-digit
        for(Map<HlaAllele,List<HlaAllele>> proteinEntry : mAlleleMap.values())
        {
            for(List<HlaAllele> alleles : proteinEntry.values())
            {
                HlaAllele matchedAllele = alleles.stream().filter(x -> x.matches(alleleStr)).findFirst().orElse(null);
                if(matchedAllele != null)
                    return matchedAllele;
            }
        }

        return null;
    }

}
