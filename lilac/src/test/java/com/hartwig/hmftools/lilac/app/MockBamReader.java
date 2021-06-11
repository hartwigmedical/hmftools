package com.hartwig.hmftools.lilac.app;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.variant.SomaticVariant;

import org.apache.commons.compress.utils.Lists;

public class MockBamReader implements BamReader
{
    public final List<Fragment> Fragments;
    public final List<Fragment> VariantFragments;
    public final Map<Indel,List<Fragment>> StopLossFragments;

    public MockBamReader()
    {
        Fragments = Lists.newArrayList();
        VariantFragments = Lists.newArrayList();
        StopLossFragments = Maps.newHashMap();
    }

    public List<Fragment> readFromBam() { return Fragments; }

    public List<Fragment> readFromBam(final SomaticVariant variant) { return VariantFragments; }

    // metrics
    public Map<Indel,List<Fragment>> getKnownStopLossFragments() { return StopLossFragments; }
    public Map<Indel,Integer> unmatchedPonIndels(int minCount) { return Maps.newHashMap(); }
    public Map<Indel,Integer> unmatchedIndels(int minCount) { return Maps.newHashMap(); }
    public int alignmentFiltered() { return 0; }

}
