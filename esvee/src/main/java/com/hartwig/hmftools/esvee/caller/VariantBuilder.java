package com.hartwig.hmftools.esvee.caller;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;

public class VariantBuilder
{
    private final HotspotCache mHotspotCache;
    private final TargetRegions mTargetRegions;
    private final StructuralVariantFactory mSvFactory;
    private final Set<String> mHotspotCandidateVcfIds;
    private final Set<String> mHardFilteredVcfIds;

    private int mHardFilteredCount;

    public VariantBuilder(final HotspotCache hotspotCache, final TargetRegions targetRegions)
    {
        mHotspotCache = hotspotCache;
        mTargetRegions = targetRegions;

        mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));
        mHardFilteredVcfIds = Sets.newHashSet();
        mHotspotCandidateVcfIds = Sets.newHashSet();
        mHardFilteredCount = 0;
    }

    public void setGenotypeOrdinals(final GenotypeIds genotypeIds)
    {
        mSvFactory.setGenotypeOrdinals(genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);
    }

    public int hardFilteredCount() { return mHardFilteredCount; }
    public int incompleteSVs() { return mSvFactory.unmatched().size(); }

    public void clearState()
    {
        mSvFactory.clear();
        mHardFilteredCount = 0;
        mHardFilteredVcfIds.clear();
        mHotspotCandidateVcfIds.clear();
    }

    public Variant checkCreateVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        if(!HumanChromosome.contains(variant.getContig()))
            return null;

        // each SV breakend can be a) not hard-filtered, b) hard-filtered but a hotspot candidate or c) neither
        // and if it's not a single then these 3 scenarios need to be considered together for the pair of breakends
        // if either are hard-filtered and not hotspot candidates, then drop them both
        boolean isSgl = StructuralVariantFactory.isSingleBreakend(variant);

        boolean hardFiltered = false;

        if(isSgl)
        {
            if(mTargetRegions.hasTargetRegions())
                hardFiltered = !mTargetRegions.inTargetRegions(variant.getContig(), variant.getStart());

            if(hardFiltered)
            {
                ++mHardFilteredCount;
                return null;
            }

            StructuralVariant sv = mSvFactory.createSingleBreakend(variant);
            return new Variant(sv, genotypeIds);
        }

        String mateId = StructuralVariantFactory.mateId(variant);

        if(mateId == null)
            return null;

        if(mHardFilteredVcfIds.contains(mateId))
        {
            // drop since the mate was hard-filtered and not a hotspot candidate, so doesn't matter what this breakend is
            mHardFilteredVcfIds.remove(mateId);
            ++mHardFilteredCount;
            return null;
        }

        boolean mateHotspotCandidate = mHotspotCandidateVcfIds.contains(mateId);

        boolean hotspotCandidate = mHotspotCache.matchesHotspotBreakend(variant.getContig(), variant.getStart());

        if(hardFiltered && !hotspotCandidate)
        {
            // clean up the first leg if it wasn't hard-filtered
            if(mateHotspotCandidate)
            {
                mHotspotCandidateVcfIds.remove(mateId);
                ++mHardFilteredCount;
            }
            else if(mSvFactory.hasUnmatchedVariant(mateId))
            {
                mSvFactory.removeUnmatchedVariant(mateId);
                ++mHardFilteredCount;
            }
            else
            {
                mHardFilteredVcfIds.add(variant.getID());
            }

            return null;
        }

        // remaining scenarios: neither leg hard-filtered, or one or both are hotspot candidates

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // check if both breakends have now been encountered
        if(currentSvCount == mSvFactory.results().size())
        {
            // this is the first breakend in the SV, so cache its state accordingly
            if(hotspotCandidate)
            {
                mHotspotCandidateVcfIds.add(variant.getID());
            }
            else if(hardFiltered)
            {
                mHardFilteredVcfIds.add(variant.getID());
            }

            return null;
        }

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return null;

        if(mateHotspotCandidate)
            mHotspotCandidateVcfIds.remove(mateId);

        if(hotspotCandidate || mateHotspotCandidate)
        {
            // check whether this SV can be rescued as a hotspot
            if(!mHotspotCache.isHotspotVariant(sv))
            {
                ++mHardFilteredCount;
                return null;
            }
        }

        return new Variant(sv, genotypeIds);
    }

    private StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

}
