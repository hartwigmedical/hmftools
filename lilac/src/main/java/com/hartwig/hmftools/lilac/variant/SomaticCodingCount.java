package com.hartwig.hmftools.lilac.variant;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.List;
import java.util.Set;

public class SomaticCodingCount
{
    public final HlaAllele Allele;

    public final double InframeIndel;
    public final double Missense;
    public final double Nonsense;
    public final double Splice;
    public final double Synonymous;
    public final double Total;

    public SomaticCodingCount(final HlaAllele allele, double inframeIndel, double missense, double nonsense, double splice, double synonymous)
    {
        Allele = allele;
        InframeIndel = inframeIndel;
        Missense = missense;
        Nonsense = nonsense;
        Splice = splice;
        Synonymous = synonymous;
        Total = InframeIndel + Missense + Nonsense + Splice + Synonymous;
    }

    public String toString()
    {
        return "SomaticCodingCount(allele=" + Allele + ", inframeIndel=" + InframeIndel + ", missense=" + Missense
                + ", nonsense=" + Nonsense + ", splice=" + Splice + ", synonymous=" + Synonymous + ")";
    }

    private static SomaticCodingCount addVariant(boolean indel, CodingEffect effect, double contribution)
    {
        SomaticCodingCount somaticCodingCount = null;


        /*
        if(indel && effect == CodingEffect.MISSENSE)
        {
            return SomaticCodingCount.copy$default(this, null, InframeIndel + contribution, 0.0, 0.0, 0.0, 0.0, 61, null);
        }
        switch(SomaticCodingCount$WhenMappings.$EnumSwitchMapping$0[effect.ordinal()])
        {
            case 1:
            {
                somaticCodingCount =
                        SomaticCodingCount.copy$default(this, null, 0.0, Missense + contribution, 0.0, 0.0, 0.0, 59, null);
                break;
            }
            case 2:
            {
                somaticCodingCount =
                        SomaticCodingCount.copy$default(this, null, 0.0, 0.0, Nonsense + contribution, 0.0, 0.0, 55, null);
                break;
            }
            case 3:
            {
                somaticCodingCount = SomaticCodingCount.copy$default(this, null, 0.0, 0.0, 0.0, Splice + contribution, 0.0, 47, null);
                break;
            }
            case 4:
            {
                somaticCodingCount =
                        SomaticCodingCount.copy$default(this, null, 0.0, 0.0, 0.0, 0.0, Synonymous + contribution, 31, null);
                break;
            }
            default:
            {
                somaticCodingCount = this;
            }
        }

         */
        return somaticCodingCount;
    }


    public static List<SomaticCodingCount> create(final List<HlaAllele> winners)
    {
        return Lists.newArrayList();

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) CollectionsKt.sorted((Iterable) winners);
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            HlaAllele hlaAllele = (HlaAllele) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            SomaticCodingCount somaticCodingCount = new SomaticCodingCount((HlaAllele) it, 0.0, 0.0, 0.0, 0.0, 0.0);
            collection.add(somaticCodingCount);
        }
        return (List) destination$iv$iv;

         */
    }

    public static List<SomaticCodingCount> addVariant(final List<SomaticCodingCount> $receiver,
            final VariantContextDecorator variant, final Set<HlaAllele> variantAlleles)
    {
        boolean bl = variant.alt().length() != variant.ref().length();
        CodingEffect codingEffect = variant.canonicalCodingEffect();
        return addVariant($receiver, bl, codingEffect, variantAlleles);
    }

    public static final List<SomaticCodingCount> addVariant(
            final List<SomaticCodingCount> codingCounts, boolean indel, final CodingEffect effect, final Set<HlaAllele> variantAlleles)
    {
        return Lists.newArrayList();

        /*
        double contribution = 1.0 / (double) variantAlleles.size();

        List result = new ArrayList();
        Iterable iterable = $receiver;
        List list = result;
        Object object = $receiver$iv;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            SomaticCodingCount it = (SomaticCodingCount) element$iv$iv;
            boolean bl = false;
            if(!(!variantAlleles.contains(it.getAllele())))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List list2 = (List) destination$iv$iv;
        list.addAll(list2);
        for(HlaAllele variantAllele : variantAlleles)
        {
            void $receiver$iv$iv2;
            Object element$iv$iv;
            Iterable $receiver$iv2 = $receiver;
            element$iv$iv = $receiver$iv2;
            Collection destination$iv$iv2 = new ArrayList();
            for(Object element$iv$iv2 : $receiver$iv$iv2)
            {
                SomaticCodingCount it = (SomaticCodingCount) element$iv$iv2;
                boolean bl = false;
                if(!Intrinsics.areEqual((Object) it.getAllele(), (Object) variantAllele))
                {
                    continue;
                }
                destination$iv$iv2.add(element$iv$iv2);
            }
            $receiver$iv2 = (List) destination$iv$iv2;
            Iterable iterable2 = $receiver$iv2;
            Comparator comparator = new Comparator<T>()
            {

                public final int compare(T a, T b)
                {
                    SomaticCodingCount it = (SomaticCodingCount) a;
                    boolean bl = false;
                    Comparable comparable = Double.valueOf(it.getTotal());
                    it = (SomaticCodingCount) b;
                    Comparable comparable2 = comparable;
                    bl = false;
                    Double d = it.getTotal();
                    return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) d);
                }
            };
            List counts = CollectionsKt.sortedWith((Iterable) iterable2, (Comparator) comparator);
            if(!(!($receiver$iv2 = (Collection) counts).isEmpty()))
            {
                continue;
            }
            result.add(((SomaticCodingCount) counts.get(0)).addVariant(indel, effect, contribution));
            result.addAll(CollectionsKt.takeLast((List) counts, (int) (counts.size() - 1)));
        }
        $receiver$iv = result;
        object = $receiver$iv;
        Comparator comparator = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                SomaticCodingCount it = (SomaticCodingCount) a;
                boolean bl = false;
                Comparable comparable = it.getAllele();
                it = (SomaticCodingCount) b;
                Comparable comparable2 = comparable;
                bl = false;
                HlaAllele hlaAllele = it.getAllele();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) hlaAllele);
            }
        };
        return CollectionsKt.sortedWith((Iterable) object, (Comparator) comparator);

         */
    }

}
