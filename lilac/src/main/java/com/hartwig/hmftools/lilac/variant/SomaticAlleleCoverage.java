package com.hartwig.hmftools.lilac.variant;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.sam.SAMRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import java.util.Collection;
import java.util.List;
import java.util.Set;

public class SomaticAlleleCoverage
{
    private final List<Integer> mVariantLoci;
    private final Set<Integer> mHetLociSansVariants;
    private final LilacConfig mConfig;
    private final Set<HlaSequenceLoci> mWinners;

    public SomaticAlleleCoverage(final LilacConfig config, final Collection<Integer> hetLoci, final LociPosition lociPosition,
            final List<? extends VariantContextDecorator> variants, final Set<HlaSequenceLoci> winners)
    {
        mConfig = config;
        mWinners = winners;
        mHetLociSansVariants = Sets.newHashSet();
        mVariantLoci = Lists.newArrayList();
        /*
        Iterable iterable = variants;
        SomaticAlleleCoverage somaticAlleleCoverage = this;
        void var7_8 = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            VariantContextDecorator variantContextDecorator = (VariantContextDecorator) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            n = lociPosition.nucelotideLoci((int) it.position());
            collection.add(n);
        }
        collection = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv = (Iterable) collection;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            int it = ((Number) element$iv$iv).intValue();
            boolean bl = false;
            if(!(it >= 0))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        collection = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv = (Iterable) collection;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            int it = ((Number) item$iv$iv).intValue();
            collection = destination$iv$iv;
            boolean bl = false;
            n = it / 3;
            collection.add(n);
        }
        collection = (List) destination$iv$iv;
        somaticAlleleCoverage.mVariantLoci = collection;
        this.mHetLociSansVariants = CollectionsKt.subtract((Iterable) hetLoci, (Iterable) this.mVariantLoci);

         */
    }

    public final List<HlaAlleleCoverage> alleleCoverage(final VariantContextDecorator variant, final SAMRecordReader reader)
    {
        return Lists.newArrayList();

        /*
        void $receiver$iv$iv;
        NucleotideFragment nucleotideFragment;
        NucleotideFragment it;
        Collection collection;
        Object item$iv$iv2;
        Iterable $receiver$iv$iv2;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull((Object) variant, (String) "variant");
        Intrinsics.checkParameterIsNotNull((Object) reader, (String) "reader");
        Iterable iterable = reader.readFromBam(variant);
        void var5_4 = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv2 : $receiver$iv$iv2)
        {
            NucleotideFragment nucleotideFragment2 = (NucleotideFragment) item$iv$iv2;
            collection = destination$iv$iv;
            boolean bl = false;
            nucleotideFragment = it.qualityFilter(this.config.getMinBaseQual());
            collection.add(nucleotideFragment);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            it = (NucleotideFragment) element$iv$iv;
            boolean bl = false;
            if(!it.isNotEmpty())
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv2 = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv2 : $receiver$iv$iv2)
        {
            it = (NucleotideFragment) item$iv$iv2;
            collection = destination$iv$iv;
            boolean bl = false;
            nucleotideFragment = it.toAminoAcidFragment();
            collection.add(nucleotideFragment);
        }
        List variantFragments = (List) destination$iv$iv;
        List<FragmentAlleles> variantFragmentAlleles =
                FragmentAlleles.Companion.create(variantFragments, (Collection<Integer>) this.hetLociSansVariants, (Collection<HlaSequenceLoci>) this.winners, (Collection<Integer>) CollectionsKt
                        .emptyList(), (Collection<HlaSequenceLoci>) CollectionsKt.emptyList());
        Iterable $receiver$iv2 = HlaAlleleCoverage.Companion.proteinCoverage(variantFragmentAlleles);
        Object object = $receiver$iv2;
        item$iv$iv2 = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                HlaAlleleCoverage it = (HlaAlleleCoverage) b;
                boolean bl = false;
                Comparable comparable = Double.valueOf(it.getTotalCoverage());
                it = (HlaAlleleCoverage) a;
                Comparable comparable2 = comparable;
                bl = false;
                Double d = it.getTotalCoverage();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) d);
            }
        };
        List coverage = CollectionsKt.sortedWith((Iterable) object, item$iv$iv2);
        $receiver$iv2 = coverage;
        object = $receiver$iv2;
        Collection destination$iv$iv2 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            HlaAlleleCoverage it2 = (HlaAlleleCoverage) element$iv$iv;
            boolean bl = false;
            if(!(it2.getTotalCoverage() == ((HlaAlleleCoverage) coverage.get(0)).getTotalCoverage()))
            {
                continue;
            }
            destination$iv$iv2.add(element$iv$iv);
        }
        return (List) destination$iv$iv2;

         */
    }

}
