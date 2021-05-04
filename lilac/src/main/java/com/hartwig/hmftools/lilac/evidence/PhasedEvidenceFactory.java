package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles;

import java.util.List;

public class PhasedEvidenceFactory
{
    private final LilacConfig mConfig;

    public PhasedEvidenceFactory(final LilacConfig config)
    {
        mConfig = config;
    }

    public final List<PhasedEvidence> evidence(final HlaContext context, final List<AminoAcidFragment> fragments)
    {
        LL_LOGGER.info("Phasing HLA-" + context.Gene + " records:");
        List<PhasedEvidence> result = this.evidence(context.ExpectedAlleles, fragments);
        if(this.mConfig.DebugPhasing)
        {
            LL_LOGGER.info("    Consolidating evidence");
        }
        for(PhasedEvidence phasedEvidence : result)
        {
            LL_LOGGER.info("    " + phasedEvidence);
        }
        return result;
    }

    public final List<PhasedEvidence> evidence(
            final ExpectedAlleles expectedAlleles, final List<AminoAcidFragment> aminoAcidAminoAcidFragments)
    {
        return Lists.newArrayList();

        /*
        Object object;
        Collection collection;
        Intrinsics.checkParameterIsNotNull((Object) expectedAlleles, (String) "expectedAlleles");
        Intrinsics.checkParameterIsNotNull(aminoAcidAminoAcidFragments, (String) "aminoAcidAminoAcidFragments");
        SequenceCount aminoAcidCounts = SequenceCount.Companion.aminoAcids(this.mConfig.getMinEvidence(), aminoAcidAminoAcidFragments);
        List<Integer> heterozygousIndices = aminoAcidCounts.heterozygousLoci();
        if(this.mConfig.getDebugPhasing())
        {
            LL_LOGGER.info("    Heterozygous Indices: " + heterozygousIndices);
        }
        ExtendEvidence heterozygousEvidence =
                new ExtendEvidence(this.mConfig, heterozygousIndices, aminoAcidAminoAcidFragments, expectedAlleles);
        Set finalisedEvidence = new LinkedHashSet();
        List unprocessedEvidence = new ArrayList();
        unprocessedEvidence.addAll((Collection) heterozygousEvidence.pairedEvidence());
        if(this.mConfig.getDebugPhasing())
        {
            LL_LOGGER.info("    Extending paired evidence");
        }
        while(!(collection = (Collection) unprocessedEvidence).isEmpty())
        {
            void parent;
            PhasedEvidence top = (PhasedEvidence) unprocessedEvidence.remove(0);
            if(this.mConfig.getDebugPhasing())
            {
                LL_LOGGER.info("    Processing top: " + top);
            }
            Object object2 = heterozygousEvidence.merge(top, SetsKt.plus((Set) finalisedEvidence, (Iterable) unprocessedEvidence));
            object = (PhasedEvidence) object2.component1();
            Set children = (Set) object2.component2();
            if(!(object2 = (Collection) children).isEmpty())
            {
                if(this.mConfig.getDebugPhasing())
                {
                    LL_LOGGER.info("    Produced child: " + parent);
                }
                finalisedEvidence.removeAll(children);
                unprocessedEvidence.removeAll(children);
                unprocessedEvidence.add(parent);
            }
            else
            {
                finalisedEvidence.add(parent);
            }
            CollectionsKt.sort((List) unprocessedEvidence);
        }
        Iterable $receiver$iv = finalisedEvidence;
        object = $receiver$iv;
        Comparator comparator = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                PhasedEvidence it = (PhasedEvidence) a;
                boolean bl = false;
                Comparable comparable = Integer.valueOf(it.getAminoAcidIndices()[0]);
                it = (PhasedEvidence) b;
                Comparable comparable2 = comparable;
                bl = false;
                Integer n = it.getAminoAcidIndices()[0];
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) n);
            }
        };
        return CollectionsKt.sortedWith((Iterable) object, (Comparator) comparator);

         */
    }
}
