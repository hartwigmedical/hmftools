package com.hartwig.hmftools.lilac.candidates;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public final class Candidates
{
    private final LilacConfig mConfig;
    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mNminoAcidSequences;

    public Candidates(final LilacConfig config, final List<HlaSequenceLoci> nucleotideSequences, final List<HlaSequenceLoci> aminoAcidSequences)
    {
        mConfig = config;
        mNucleotideSequences = nucleotideSequences;
        mNminoAcidSequences = aminoAcidSequences;
    }

    public final List<HlaAllele> unphasedCandidates(final HlaContext context, final List<AminoAcidFragment> fragments)
    {
        return Lists.newArrayList();

        /*
        void $receiver$iv$iv;
        void $receiver$iv;
        Iterable $receiver$iv$iv2;
        void $receiver$iv2;
        Iterable $receiver$iv$iv3;
        HlaAllele hlaAllele;
        Collection collection;
        Object item$iv$iv2;
        void $receiver$iv$iv4;
        Object element$iv$iv2;
        void $receiver$iv$iv5;
        Iterable $receiver$iv3;
        Intrinsics.checkParameterIsNotNull((Object) context, (String) "context");
        Intrinsics.checkParameterIsNotNull(fragments, (String) "fragments");
        String gene = context.getGene();
        Set<Integer> aminoAcidBoundary = context.getAminoAcidBoundaries();
        logger.info("Determining un-phased candidate set for gene HLA-" + gene);
        SequenceCount aminoAcidCounts = SequenceCount.Companion.aminoAcids(mConfig.getMinEvidence(), fragments);
        Iterable iterable = $receiver$iv3 = (Iterable) mNminoAcidSequences;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv2 : $receiver$iv$iv5)
        {
            HlaSequenceLoci it = (HlaSequenceLoci) element$iv$iv2;
            boolean bl = false;
            if(!Intrinsics.areEqual((Object) it.getAllele().getGene(), (Object) gene))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv2);
        }
        List geneCandidates = (List) destination$iv$iv;
        logger.info("    " + geneCandidates.size() + " candidates before filtering");
        com.hartwig.hmftools.lilackt.candidates.AminoAcidFiltering aminoAcidFilter = new AminoAcidFiltering(aminoAcidBoundary);
        List<HlaSequenceLoci> aminoAcidCandidates = aminoAcidFilter.aminoAcidCandidates(geneCandidates, aminoAcidCounts);
        Iterable $receiver$iv4 = aminoAcidCandidates;
        element$iv$iv2 = $receiver$iv4;
        Iterable destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv4, (int) 10));
        for(Object item$iv$iv2 : $receiver$iv$iv4)
        {
            void it;
            HlaSequenceLoci $i$f$filter = (HlaSequenceLoci) item$iv$iv2;
            collection = destination$iv$iv2;
            boolean bl = false;
            hlaAllele = it.getAllele();
            collection.add(hlaAllele);
        }
        Set aminoAcidCandidateAlleles = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv2));
        Iterable $receiver$iv5 = aminoAcidCandidateAlleles;
        destination$iv$iv2 = $receiver$iv5;
        Iterable destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv5, (int) 10));
        for(Object item$iv$iv3 : $receiver$iv$iv3)
        {
            void it;
            HlaAllele bl = (HlaAllele) item$iv$iv3;
            collection = destination$iv$iv3;
            boolean bl2 = false;
            hlaAllele = it.asFourDigit();
            collection.add(hlaAllele);
        }
        Set aminoAcidSpecificAllelesCandidate = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv3));
        if(aminoAcidSpecificAllelesCandidate.isEmpty())
        {
            logger.warn("    0 candidates after amino acid filtering - reverting to all gene candidates");
            $receiver$iv$iv3 = $receiver$iv5 = (Iterable) geneCandidates;
            destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv5, (int) 10));
            for(Object item$iv$iv3 : $receiver$iv$iv3)
            {
                HlaSequenceLoci it = (HlaSequenceLoci) item$iv$iv3;
                collection = destination$iv$iv3;
                boolean bl = false;
                hlaAllele = it.getAllele();
                collection.add(hlaAllele);
            }
            return (List) destination$iv$iv3;
        }
        logger.info("    " + aminoAcidCandidates.size() + " candidates after amino acid filtering");
        com.hartwig.hmftools.lilackt.candidates.NucleotideFiltering
                nucleotideFiltering = new NucleotideFiltering(mConfig.getMinEvidence(), aminoAcidBoundary);
        destination$iv$iv3 = mNucleotideSequences;
        item$iv$iv2 = $receiver$iv2;
        Collection destination$iv$iv4 = new ArrayList();
        for(Object element$iv$iv3 : $receiver$iv$iv2)
        {
            HlaSequenceLoci it = (HlaSequenceLoci) element$iv$iv3;
            boolean bl = false;
            if(!aminoAcidSpecificAllelesCandidate.contains(it.getAllele().asFourDigit()))
            {
                continue;
            }
            destination$iv$iv4.add(element$iv$iv3);
        }
        List nucleotideCandidatesAfterAminoAcidFiltering = (List) destination$iv$iv4;
        $receiver$iv$iv2 =
                nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(nucleotideCandidatesAfterAminoAcidFiltering, fragments);
        destination$iv$iv4 = $receiver$iv;
        Collection destination$iv$iv5 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv4 : $receiver$iv$iv)
        {
            void it;
            HlaSequenceLoci bl = (HlaSequenceLoci) item$iv$iv4;
            collection = destination$iv$iv5;
            boolean bl3 = false;
            hlaAllele = it.getAllele().asFourDigit();
            collection.add(hlaAllele);
        }
        List nucleotideSpecificAllelesCandidate = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv5));
        if(nucleotideSpecificAllelesCandidate.isEmpty())
        {
            logger.warn("    0 candidates after exon boundary filtering - reverting to amino acid candidates");
            return CollectionsKt.toList((Iterable) aminoAcidCandidateAlleles);
        }
        logger.info("    " + nucleotideSpecificAllelesCandidate.size() + " candidates after exon boundary filtering");
        return nucleotideSpecificAllelesCandidate;

         */
    }

    public final List<HlaAllele> phasedCandidates(
            final HlaContext context, final List<HlaAllele> unphasedCandidateAlleles, final List<PhasedEvidence> phasedEvidence)
    {
        return Lists.newArrayList();

        /*
        Collection<HlaAllele> collection;
        Iterable $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable $receiver$iv$iv2;
        Iterable $receiver$iv2;
        Intrinsics.checkParameterIsNotNull((Object) context, (String) "context");
        Intrinsics.checkParameterIsNotNull(unphasedCandidateAlleles, (String) "unphasedCandidateAlleles");
        Intrinsics.checkParameterIsNotNull(phasedEvidence, (String) "phasedEvidence");
        String gene = context.getGene();
        logger.info("Determining phased candidate set for gene HLA-" + gene);
        Iterable iterable = $receiver$iv2 = (Iterable) mNminoAcidSequences;
        Collection destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            HlaSequenceLoci it = (HlaSequenceLoci) element$iv$iv;
            boolean bl = false;
            if(!unphasedCandidateAlleles.contains(it.getAllele().asFourDigit()))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List unphasedCandidates = (List) destination$iv$iv;
        List<HlaSequenceLoci> phasedCandidates = filterCandidates(unphasedCandidates, phasedEvidence);
        $receiver$iv$iv2 = phasedCandidates;
        Comparable<StringBuilder> comparable =
                new StringBuilder().append("    ").append(phasedCandidates.size()).append(" candidates after phasing: ");
        Object object = logger;
        destination$iv$iv = $receiver$iv;
        Collection destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            HlaSequenceLoci bl = (HlaSequenceLoci) item$iv$iv;
            collection = destination$iv$iv2;
            boolean bl2 = false;
            HlaAllele hlaAllele = it.getAllele();
            collection.add(hlaAllele);
        }
        collection = (List) destination$iv$iv2;
        object.info(comparable.append(CollectionsKt.joinToString$default((Iterable) collection, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null))
                .toString());
        $receiver$iv$iv = $receiver$iv = (Iterable) phasedCandidates;
        destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            HlaSequenceLoci it = (HlaSequenceLoci) item$iv$iv;
            object = destination$iv$iv2;
            boolean bl = false;
            comparable = it.getAllele();
            object.add(comparable);
        }
        return (List) destination$iv$iv2;

         */
    }

    private List<HlaSequenceLoci> filterCandidates(final List<HlaSequenceLoci> initialCandidates, final List<PhasedEvidence> evidence)
    {
        /*
        var candidates = initialCandidates
        for (i in evidence.indices) {
            val newEvidence = evidence[i]
            candidates = candidates.filter { it.consistentWithAny(newEvidence.evidence.keys, *newEvidence.aminoAcidIndices) }
        }

        return candidates
         */

        List<HlaSequenceLoci> candidates = Lists.newArrayList();
        candidates.addAll(initialCandidates);

        for(int i = 0; i < evidence.size(); ++i)
        {
            PhasedEvidence newEvidence = evidence.get(i);

            candidates = candidates.stream()
                .filter(x -> x.consistentWithAny(newEvidence.getEvidence().keySet(), newEvidence.getAminoAcidIndices()))
                .collect(Collectors.toList());
        }

        return candidates;
    }

}
