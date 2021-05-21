package com.hartwig.hmftools.lilac.variant;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragment;
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
        List<NucleotideFragment> fragments = reader.readFromBam(variant);

        List<AminoAcidFragment> variantFragments = fragments.stream()
                .map(x -> x.qualityFilter(mConfig.MinBaseQual))
                .filter(x -> x.isNotEmpty())
                .map(x -> x.toAminoAcidFragment())
                .collect(Collectors.toList());

        List<FragmentAlleles> variantFragmentAlleles = FragmentAlleles.createFragmentAlleles(
                variantFragments, mHetLociSansVariants, mWinners, Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());

        List<HlaAlleleCoverage> coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles);

        Collections.sort(coverage, new HlaAlleleCoverage.TotalCoverageSorter());

        return coverage.stream().filter(x -> x.TotalCoverage == coverage.get(0).TotalCoverage).collect(Collectors.toList());
    }
}
