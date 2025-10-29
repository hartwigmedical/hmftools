package com.hartwig.hmftools.pavereverse.parse;

import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSpecification;
import com.hartwig.hmftools.pavereverse.protein.Deletion;
import com.hartwig.hmftools.pavereverse.protein.DeletionInsertion;
import com.hartwig.hmftools.pavereverse.protein.Duplication;
import com.hartwig.hmftools.pavereverse.protein.Frameshift;
import com.hartwig.hmftools.pavereverse.protein.Insertion;
import com.hartwig.hmftools.pavereverse.protein.ProteinVariant;
import com.hartwig.hmftools.pavereverse.protein.SingleAminoAcidVariant;
import com.hartwig.hmftools.pavereverse.protein.StartLost;
import com.hartwig.hmftools.pavereverse.protein.StopGained;

class ProteinVariantFactory
{
    private final GeneData mGeneData;
    private final AminoAcidRange mRefRange;
    private final Set<ProteinTranscript> mTranscripts;

    ProteinVariantFactory(GeneData geneData, AminoAcidRange refRange, Set<ProteinTranscript> transcript)
    {
        mGeneData = geneData;
        mRefRange = refRange;
        mTranscripts = transcript;
    }

    Set<ProteinVariant> buildSingleAminoAcidVariants(AminoAcidSpecification altSpec)
    {
        return mTranscripts.stream()
                .map(t -> new SingleAminoAcidVariant(mGeneData, t.mTranscriptData, t.mAminoAcids, altSpec))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildDuplication()
    {
        return mTranscripts.stream()
                .map(t -> new Duplication(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildFrameshift()
    {
        return mTranscripts.stream()
                .map(t -> new Frameshift(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildStopGained()
    {
        return mTranscripts.stream()
                .map(t -> new StopGained(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildStartLost()
    {
        return mTranscripts.stream()
                .map(t -> new StartLost(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildDeletion()
    {
        return mTranscripts.stream()
                .map(t -> new Deletion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildDeletionInsertion(String altSequence)
    {
        AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(altSequence);
        return mTranscripts.stream()
                .map(t -> new DeletionInsertion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange, altAminoAcids))
                .collect(Collectors.toSet());
    }

    Set<ProteinVariant> buildInsertion(String inserted)
    {
        AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(inserted);
        return mTranscripts.stream()
                .map(t -> new Insertion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange, altAminoAcids))
                .collect(Collectors.toSet());
    }
}
