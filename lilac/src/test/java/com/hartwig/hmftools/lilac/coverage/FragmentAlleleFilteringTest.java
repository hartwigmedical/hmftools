package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.filterUnsupportedWildcardFragments;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findUnsupportedWildcards;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper.findWildcardAlleles;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createFragment;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class FragmentAlleleFilteringTest
{
    @Test
    public void testWildcardFiltering()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaAllele allele2 = HlaAllele.fromString("A*01:02");

        HlaAllele wcAllele3 = HlaAllele.fromString("A*01:56");
        wcAllele3.setHasWildcard(true);

        HlaAllele wcAllele4 = HlaAllele.fromString("A*01:77");
        wcAllele4.setHasWildcard(true);

        HlaAllele wcAllele5 = HlaAllele.fromString("A*01:102");
        wcAllele5.setHasWildcard(true);

        List<FragmentAlleles> fragAlleles = Lists.newArrayList();

        // the first wc allele will have sufficient fragments, the second not due to lack of frags, the third not due to non-unique frags
        fragAlleles.add(new FragmentAlleles(createFragment("01"), Lists.newArrayList(wcAllele3), Lists.newArrayList()));
        fragAlleles.add(new FragmentAlleles(createFragment("02"), Lists.newArrayList(wcAllele3, wcAllele4, wcAllele5), Lists.newArrayList()));

        // non-unique doesn't count
        fragAlleles.add(new FragmentAlleles(createFragment("03"), Lists.newArrayList(allele1), Lists.newArrayList(wcAllele5)));
        fragAlleles.add(new FragmentAlleles(createFragment("04"), Lists.newArrayList(allele2), Lists.newArrayList(wcAllele5)));

        // wild-only doesn't count
        fragAlleles.add(new FragmentAlleles(createFragment("05"), Lists.newArrayList(), Lists.newArrayList(wcAllele3, wcAllele4, wcAllele5)));

        // will be filtered out once unsupported WCs are done
        fragAlleles.add(new FragmentAlleles(createFragment("06"), Lists.newArrayList(), Lists.newArrayList(wcAllele4)));
        fragAlleles.add(new FragmentAlleles(createFragment("07"), Lists.newArrayList(), Lists.newArrayList(wcAllele5)));
        fragAlleles.add(new FragmentAlleles(createFragment("08"), Lists.newArrayList(), Lists.newArrayList(wcAllele4, wcAllele5)));

        Set<HlaAllele> wildcardAlleles = findWildcardAlleles(fragAlleles);
        assertEquals(3, wildcardAlleles.size());

        List<HlaAllele> unsupported = findUnsupportedWildcards(fragAlleles, wildcardAlleles);
        assertEquals(2, unsupported.size());
        assertTrue(unsupported.contains(wcAllele4));
        assertTrue(unsupported.contains(wcAllele5));

        filterUnsupportedWildcardFragments(fragAlleles, unsupported);
        assertEquals(5, fragAlleles.size());
        assertFalse(fragAlleles.stream().anyMatch(x -> x.getFragment().id().equals("06")));
        assertFalse(fragAlleles.stream().anyMatch(x -> x.getFragment().id().equals("07")));
        assertFalse(fragAlleles.stream().anyMatch(x -> x.getFragment().id().equals("08")));
    }
}
