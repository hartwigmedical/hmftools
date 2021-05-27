package com.hartwig.hmftools.lilac.variant;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class SomaticAlleleCoverage
{
    private final List<Integer> mVariantLoci;
    private final List<Integer> mHetLociSansVariants;
    private final LilacConfig mConfig;
    private final List<HlaSequenceLoci> mWinners;

    public SomaticAlleleCoverage(final LilacConfig config, final List<Integer> hetLoci, final LociPosition lociPosition,
            final List<VariantContextDecorator> variants, final List<HlaSequenceLoci> winners)
    {
        mConfig = config;
        mWinners = winners;

        mVariantLoci = Lists.newArrayList();

        variants.stream()
                .map(x -> lociPosition.nucelotideLoci((int)x.position()))
                .filter(x -> x >= 0)
                .mapToInt(x -> x / 3)
                .forEach(x -> mVariantLoci.add(x));

        mHetLociSansVariants = hetLoci.stream().filter(x -> !mVariantLoci.contains(x)).collect(Collectors.toList());
    }

    public final List<HlaAlleleCoverage> alleleCoverage(final VariantContextDecorator variant, final SAMRecordReader reader)
    {
        List<Fragment> fragments = reader.readFromBam(variant);

        fragments.forEach(x -> x.qualityFilter(mConfig.MinBaseQual));
        fragments.forEach(x -> x.buildAminoAcids());

        List<Fragment> variantFragments = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            fragment.qualityFilter(mConfig.MinBaseQual);
            fragment.buildAminoAcids();

            if(fragment.hasNucleotides())
                variantFragments.add(fragment);
        }

        List<FragmentAlleles> variantFragmentAlleles = FragmentAlleles.createFragmentAlleles(
                variantFragments, mHetLociSansVariants, mWinners, Lists.newArrayList(), Maps.newHashMap(),
                Lists.newArrayList(), Lists.newArrayList());

        List<HlaAlleleCoverage> coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles);

        Collections.sort(coverage, new HlaAlleleCoverage.TotalCoverageSorter());

        return coverage.stream().filter(x -> x.TotalCoverage == coverage.get(0).TotalCoverage).collect(Collectors.toList());
    }
}
