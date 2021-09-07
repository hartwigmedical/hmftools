package com.hartwig.hmftools.lilac.read;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.variant.SomaticVariant;

public interface BamReader
{
    List<Fragment> readFromBam();

    List<Fragment> readFromBam(final SomaticVariant variant);

    // metrics
    Map<Indel,List<Fragment>> getKnownStopLossFragments();
    Map<Indel,Integer> unmatchedPonIndels(int minCount);
    Map<Indel,Integer> unmatchedIndels(int minCount);
    int alignmentFiltered();
}
