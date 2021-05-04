package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class PhasedEvidenceValidation
{
    public static void validateExpected(final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        /*
        void $receiver$iv$iv;
        Intrinsics.checkParameterIsNotNull((Object) gene, (String) "gene");
        Intrinsics.checkParameterIsNotNull(evidence, (String) "evidence");
        Intrinsics.checkParameterIsNotNull(candidates, (String) "candidates");
        Iterable $receiver$iv = candidates;
        Iterable iterable = $receiver$iv;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            HlaSequenceLoci it = (HlaSequenceLoci) element$iv$iv;
            boolean bl = false;
            if(!Intrinsics.areEqual((Object) it.getAllele().getGene(), (Object) gene))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List expectedSequences = (List) destination$iv$iv;
        for(HlaSequenceLoci sequence2 : expectedSequences)
        {
            for(PhasedEvidence phasedEvidence : evidence)
            {
                if(sequence2.consistentWith(phasedEvidence))
                {
                    continue;
                }
                logger.warn("Expected allele " + sequence2.getAllele() + " filtered by " + phasedEvidence);
            }
        }
        
         */
    }

    public static void validateAgainstFinalCandidates(
            final String gene, final List<PhasedEvidence> evidence, final List<HlaSequenceLoci> candidates)
    {
        for(PhasedEvidence inconsistentEvidence2 : unmatchedEvidence(evidence, candidates))
        {
            LL_LOGGER.warn("HLA-" + gene + " phased evidence not found in candidates: " + inconsistentEvidence2);
        }
    }

    private static List<PhasedEvidence> unmatchedEvidence(List<PhasedEvidence> evidence, List<HlaSequenceLoci> candidates)
    {
        return Lists.newArrayList();

        /*
        PhasedEvidence it;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) evidence;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            PhasedEvidence phasedEvidence = (PhasedEvidence) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            PhasedEvidence phasedEvidence2 = it.inconsistentEvidence((Collection<HlaSequenceLoci>) candidates);
            collection.add(phasedEvidence2);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (PhasedEvidence) element$iv$iv;
            boolean bl = false;
            Map<String, Integer> map = it.getEvidence();
            if(!(!map.isEmpty()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        return (List) destination$iv$iv;

         */
    }
}
