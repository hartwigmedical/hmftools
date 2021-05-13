package com.hartwig.hmftools.lilac.hla;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

public class HlaAlleleCache
{
    // ensure the same HlaAlle is not created twice, allowing for object matching instead of string matching
    private final List<HlaAllele> mGroups;
    private final List<HlaAllele> mFourDigitList;
    private final List<HlaAllele> mAlleles;

    public HlaAlleleCache()
    {
        mGroups = Lists.newArrayList();
        mFourDigitList = Lists.newArrayList();
        mAlleles = Lists.newArrayList();
    }

    public int alleleCount() { return mAlleles.size(); }
    public int fourDigitCount() { return mFourDigitList.size(); }
    public int groupCount() { return mGroups.size(); }

    public HlaAllele request(final String alleleStr)
    {
        // first check if the allele has already been created
        HlaAllele match = mAlleles.stream().filter(x -> x.matches(alleleStr)).findFirst().orElse(null);
        if(match != null)
            return match;

        HlaAllele tmpAllele = HlaAllele.fromString(alleleStr);

        // use existing group and protein/4-digit versions if they exist already
        HlaAllele group = requestGroup(tmpAllele.asAlleleGroup());
        HlaAllele protein = requestFourDigit(tmpAllele.asFourDigit().toString());

        HlaAllele newAllele = new HlaAllele(
                tmpAllele.Gene, tmpAllele.AlleleGroup, tmpAllele.Protein, tmpAllele.Synonymous, tmpAllele.SynonymousNonCoding,
                protein, group);

        mAlleles.add(newAllele);
        return newAllele;
    }

    public HlaAllele requestFourDigit(final String alleleStr)
    {
        HlaAllele tmpAllele = HlaAllele.fromString(alleleStr).asFourDigit();

        HlaAllele match = mFourDigitList.stream().filter(x -> x.matches(tmpAllele)).findFirst().orElse(null);
        if(match != null)
            return match;

        HlaAllele group = requestGroup(tmpAllele.asAlleleGroup());

        HlaAllele newAllele = new HlaAllele(
                tmpAllele.Gene, tmpAllele.AlleleGroup, tmpAllele.Protein, tmpAllele.Synonymous, tmpAllele.SynonymousNonCoding,
                null, group);

        mFourDigitList.add(newAllele);
        return newAllele;
    }

    public void rebuildProteinAlleles(final List<HlaAllele> list)
    {
        List<HlaAllele> cachedList = list.stream().map(x -> requestFourDigit(x.asFourDigit().toString())).collect(Collectors.toList());
        list.clear();
        list.addAll(cachedList);
    }


    /*
    public HlaAllele request(final HlaAllele allele)
    {
        HlaAllele match = mAlleles.stream().filter(x -> x.matches(allele)).findFirst().orElse(null);
        if(match != null)
            return match;

        mAlleles.add(allele);
        return allele;
    }


    public HlaAllele requestFourDigit(final String alleleStr)
    {
        return requestFourDigit(HlaAllele.fromString(alleleStr));
    }

    public HlaAllele requestGroup(final String alleleStr)
    {
        return requestGroup(HlaAllele.fromString(alleleStr));
    }
    */

    public HlaAllele requestGroup(final HlaAllele requested)
    {
        // ensure in group format
        HlaAllele allele = requested.asAlleleGroup().matches(requested) ? requested : requested.asAlleleGroup();

        HlaAllele match = mGroups.stream().filter(x -> x.matches(allele)).findFirst().orElse(null);
        if(match != null)
            return match;

        mGroups.add(allele);
        return allele;
    }

    @VisibleForTesting
    public List<HlaAllele> getAlleles() { return mAlleles; }

    @VisibleForTesting
    public List<HlaAllele> getFourDigits() { return mFourDigitList; }

    @VisibleForTesting
    public List<HlaAllele> getGroups() { return mGroups; }


}
