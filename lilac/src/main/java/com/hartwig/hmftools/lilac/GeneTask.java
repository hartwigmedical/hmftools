package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.Candidates;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class GeneTask implements Callable<Void>
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;
    private final AminoAcidFragmentPipeline mAminoAcidPipeline;
    private final Candidates mCandidateFactory;

    // gene specific information
    private final HlaContext mHlaContext;

    private final List<Fragment> mCandidateFrags;
    private final List<HlaAllele> mCandidatesAlleles;
    private final List<PhasedEvidence> mPhasedEvidence;

    public GeneTask(
            final LilacConfig config, final ReferenceData referenceData, final AminoAcidFragmentPipeline aminoAcidPipeline,
            final Candidates candidateFactory, final HlaContext hlaContext)
    {
        mConfig = config;
        mRefData = referenceData;
        mAminoAcidPipeline = aminoAcidPipeline;
        mCandidateFactory = candidateFactory;

        mHlaContext = hlaContext;

        // extract fragments relevant for this gene, copying fragment data in the process to avoid concurrent access as they are filtered
        mCandidateFrags = Lists.newArrayList(mAminoAcidPipeline.referencePhasingFragments(mHlaContext));
        mCandidatesAlleles = Lists.newArrayList();
        mPhasedEvidence = Lists.newArrayList();
    }

    public List<PhasedEvidence> phasedEvidence() { return mPhasedEvidence; }

    @Override
    public Void call()
    {
        // determine un-phased Candidates
        List<HlaAllele> unphasedCandidates = mCandidateFactory.unphasedCandidates(mHlaContext, mCandidateFrags, mRefData.CommonAlleles);

        // determine phasing of amino acids
        PhasedEvidenceFactory phasedEvidenceFactory = new PhasedEvidenceFactory(mConfig);
        mPhasedEvidence.addAll(phasedEvidenceFactory.evidence(mHlaContext, mCandidateFrags));

        // validate phasing against expected sequences
        if(!mConfig.ActualAlleles.isEmpty() && LL_LOGGER.isDebugEnabled())
        {
            List<HlaSequenceLoci> actualSequences = mRefData.AminoAcidSequences.stream()
                    .filter(x -> mConfig.ActualAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

            PhasedEvidence.logInconsistentEvidence(mHlaContext.Gene, mPhasedEvidence, actualSequences);
        }

        // gather all phased candidates
        mCandidatesAlleles.addAll(mCandidateFactory.phasedCandidates(mHlaContext, Sets.newHashSet(unphasedCandidates), mPhasedEvidence));

        return null;
    }

    public void addPhasedCandidates(final List<HlaAllele> allAlleles)
    {
        if(mCandidatesAlleles.isEmpty())
            return;

        if(mConfig.MaxEliminationCandidates == 0 || mCandidatesAlleles.size() <= mConfig.MaxEliminationCandidates)
        {
            allAlleles.addAll(mCandidatesAlleles);
            return;
        }

        final HlaGene gene = mCandidatesAlleles.get(0).Gene;

        mRefData.getAlleleFrequencies().getAlleleFrequencies().keySet().stream()
                .filter(x -> x.Gene == gene).forEach(allAlleles::add);
    }
}
