package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.ExpectedAlleles;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

import org.apache.commons.math3.util.Pair;

public class PhasedEvidenceFactory
{
    private final LilacConfig mConfig;
    private final boolean mDebugPhasing;

    public PhasedEvidenceFactory(final LilacConfig config)
    {
        mConfig = config;
        mDebugPhasing = mConfig.DebugPhasing;
    }

    public List<PhasedEvidence> evidence(final HlaContext context, final List<Fragment> fragments)
    {
        if(mDebugPhasing)
        {
            LL_LOGGER.debug("phasing {} records:", context.geneName());
        }

        List<PhasedEvidence> result = evidence(context.ExpectedAlleles, fragments);

        if(mDebugPhasing)
        {
            LL_LOGGER.debug("  consolidating evidence");
            for(PhasedEvidence phasedEvidence : result)
            {
                LL_LOGGER.debug("  " + phasedEvidence);
            }
        }
        return result;
    }

    public List<PhasedEvidence> evidence(final ExpectedAlleles expectedAlleles, final List<Fragment> fragments)
    {
        SequenceCount aminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidenceSupport, mConfig.MinEvidenceFactor, fragments);

        List<Integer> heterozygousIndices = Lists.newArrayList(aminoAcidCounts.heterozygousLoci());

        if(mDebugPhasing)
        {
            LL_LOGGER.debug("  heterozygous Indices: {}", heterozygousIndices);
        }

        ExtendEvidence heterozygousEvidence = new ExtendEvidence(mConfig, heterozygousIndices, fragments, expectedAlleles);

        List<PhasedEvidence> finalisedEvidence = Lists.newArrayList();
        List<PhasedEvidence> unprocessedEvidence = Lists.newArrayList();
        unprocessedEvidence.addAll(heterozygousEvidence.pairedEvidence());

        if(mDebugPhasing)
        {
            LL_LOGGER.debug("  extending paired evidence");
        }

        while(!unprocessedEvidence.isEmpty())
        {
            PhasedEvidence top = unprocessedEvidence.remove(0);

            if(mDebugPhasing)
            {
                LL_LOGGER.debug("  Processing top: {}", top);
            }

            Set<PhasedEvidence> others = Sets.newHashSet();
            others.addAll(finalisedEvidence);
            others.addAll(unprocessedEvidence);
            Pair<PhasedEvidence, Set<PhasedEvidence>> pair = heterozygousEvidence.merge(top, others);
            PhasedEvidence parent = pair.getFirst();
            Set<PhasedEvidence> children = pair.getSecond();

            if(!children.isEmpty())
            {
                if(mDebugPhasing)
                {
                    LL_LOGGER.debug("  Produced child: {}", pair.getFirst());
                }

                finalisedEvidence.removeAll(children);
                unprocessedEvidence.removeAll(children);
                unprocessedEvidence.add(parent);
            }
            else
            {
                finalisedEvidence.add(parent);
            }

            Collections.sort(unprocessedEvidence);
        }

        Collections.sort(finalisedEvidence, new PhasedEvidenceSorter());
        return finalisedEvidence;
    }

    private static class PhasedEvidenceSorter implements Comparator<PhasedEvidence>
    {
        public int compare(final PhasedEvidence first, final PhasedEvidence second)
        {
            int firstAA = first.getAminoAcidLoci().get(0);
            int secondAA = second.getAminoAcidLoci().get(0);
            if(firstAA != secondAA)
            {
                return firstAA > secondAA ? 1 : -1;
            }

            return 0;
        }
    }

}
